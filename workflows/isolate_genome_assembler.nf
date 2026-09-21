include { READ_QC       } from '../subworkflows/local/read_qc'
include { CONTAMINATION } from '../subworkflows/local/contamination'
include { ASSEMBLY      } from '../subworkflows/local/assembly'
include { FINISHING     } from '../subworkflows/local/finishing'
include { CHECKS        } from '../subworkflows/local/checks'
include { REMOVE_HUMAN  } from '../modules/local/remove_human'
include { QC_GATES      } from '../modules/local/qc_gates'

// The thresholds YAML as JSON, read once here so that a malformed file stops the run at
// start-up rather than in the last task of every sample.
def qcThresholdsJson(path) {
    def thresholds = new org.yaml.snakeyaml.Yaml().load(file(path, checkIfExists: true).text)
    if (!(thresholds instanceof Map) || thresholds.values().any { !(it instanceof Map) }) {
        error "Invalid --qc_thresholds '${path}': expected a mapping of check name to {direction, warn, fail}"
    }
    groovy.json.JsonOutput.toJson(thresholds)
}

workflow ISOLATE_GENOME_ASSEMBLER {
    take:
    ch_samples // tuple: [meta map, one or more PacBio HiFi read files]

    main:
    READ_QC(ch_samples)
    CONTAMINATION(READ_QC.out.reads)

    // Always runs so the DAG stays static and resume-safe; the decision file, not the
    // graph, says whether any reads are actually dropped.
    REMOVE_HUMAN(
        READ_QC.out.reads
            .join(CONTAMINATION.out.human_read_ids)
            .join(CONTAMINATION.out.remove_human)
    )

    ASSEMBLY(REMOVE_HUMAN.out.reads, READ_QC.out.summary)
    FINISHING(ASSEMBLY.out.assembly, REMOVE_HUMAN.out.reads)
    CHECKS(
        FINISHING.out.assembly,
        FINISHING.out.contigs,
        REMOVE_HUMAN.out.reads,
        READ_QC.out.summary,
        ASSEMBLY.out.autocycler_dir,
        ASSEMBLY.out.consensus_gfa,
    )

    // Every measurement for a sample, gathered into parallel lists of qc_gates.py option
    // names and files. groupTuple waits for all of a sample's checks, which is the point.
    ch_metrics = channel.empty()
        .mix(
            READ_QC.out.summary.map { meta, f -> [meta, 'read-qc', f] },
            CONTAMINATION.out.summary.map { meta, f -> [meta, 'contamination', f] },
            ASSEMBLY.out.assembly_source.map { meta, f -> [meta, 'assembly-source', f] },
            FINISHING.out.contigs.map { meta, f -> [meta, 'contigs', f] },
            FINISHING.out.plasmid_audit.map { meta, f -> [meta, 'plasmid-audit', f] },
            CHECKS.out.metrics,
        )
        .groupTuple()

    QC_GATES(ch_metrics, params.input ? qcThresholdsJson(params.qc_thresholds) : '{}')

    emit:
    reads = REMOVE_HUMAN.out.reads
    assembly = FINISHING.out.assembly
    contigs = FINISHING.out.contigs
    plasmid_audit = FINISHING.out.plasmid_audit
    assembly_source = ASSEMBLY.out.assembly_source
    read_qc = READ_QC.out.summary
    contamination = CONTAMINATION.out.summary
    qc = QC_GATES.out.qc
}

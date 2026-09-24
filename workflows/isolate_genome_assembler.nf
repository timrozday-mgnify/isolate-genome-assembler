include { READ_QC       } from '../subworkflows/local/read_qc'
include { CONTAMINATION } from '../subworkflows/local/contamination'
include { ASSEMBLY      } from '../subworkflows/local/assembly'
include { FINISHING     } from '../subworkflows/local/finishing'
include { SCORING       } from '../subworkflows/local/scoring'
include { CHECKS        } from '../subworkflows/local/checks'
include { REMOVE_HUMAN  } from '../modules/local/remove_human'
include { QC_GATES      } from '../modules/local/qc_gates'
include { COLLECT_METRICS } from '../modules/local/collect_metrics'
include { QUARTO_REPORT } from '../modules/local/quarto_report'

// The thresholds YAML as JSON, read once here so that a malformed file stops the run at
// start-up rather than in the last task of every sample.
def qcThresholdsJson(path) {
    def thresholds = new org.yaml.snakeyaml.Yaml().load(file(path, checkIfExists: true).text)
    if (!(thresholds instanceof Map) || thresholds.values().any { !(it instanceof Map) }) {
        error "Invalid --qc_thresholds '${path}': expected a mapping of check name to {direction, warn, fail}"
    }
    groovy.json.JsonOutput.toJson(thresholds)
}

// Params as JSON values: paths and GStrings become plain strings.
def jsonValue(value) {
    (value == null || value instanceof Number || value instanceof Boolean) ? value : value.toString()
}

// What the report's run overview shows: the pipeline, this run, and every param.
def runInfoJson() {
    groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        pipeline: workflow.manifest.name,
        pipeline_version: workflow.manifest.version,
        commit: workflow.commitId ?: '',
        revision: workflow.revision ?: '',
        nextflow_version: nextflow.version.toString(),
        run_name: workflow.runName,
        start: workflow.start.toString(),
        profile: workflow.profile,
        container_engine: workflow.containerEngine ?: '',
        params: params.collectEntries { key, value -> [key, jsonValue(value)] }.sort(),
    ]))
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

    // Every candidate assembly scored against the others -- each assembler's full-read
    // assembly and the consensus -- and the winner picked, as --assembly_selection says.
    SCORING(
        ASSEMBLY.out.full_assemblies.mix(
            ASSEMBLY.out.consensus.map { meta, assembly -> [meta, 'autocycler', assembly] }
        ),
        REMOVE_HUMAN.out.reads,
        READ_QC.out.summary,
        ASSEMBLY.out.selection_in,
    )

    FINISHING(SCORING.out.assembly, REMOVE_HUMAN.out.reads)
    CHECKS(
        FINISHING.out.assembly,
        FINISHING.out.contigs,
        REMOVE_HUMAN.out.reads,
        ASSEMBLY.out.autocycler_dir,
        ASSEMBLY.out.consensus_gfa,
        SCORING.out.meryl_db,
    )

    // Every measurement for a sample, gathered into parallel lists of qc_gates.py option
    // names and files. groupTuple waits for all of a sample's checks, which is the point.
    ch_measurements = channel.empty()
        .mix(
            READ_QC.out.summary.map { meta, f -> [meta, 'read-qc', f] },
            CONTAMINATION.out.summary.map { meta, f -> [meta, 'contamination', f] },
            SCORING.out.assembly_source.map { meta, f -> [meta, 'assembly-source', f] },
            FINISHING.out.contigs.map { meta, f -> [meta, 'contigs', f] },
            FINISHING.out.plasmid_audit.map { meta, f -> [meta, 'plasmid-audit', f] },
            CHECKS.out.metrics,
        )

    QC_GATES(ch_measurements.groupTuple(), params.input ? qcThresholdsJson(params.qc_thresholds) : '{}')

    // --- Stage 8: report ---
    // Every file the report reads, as [sample id, kind, file]; kinds are explained in
    // bin/collect_metrics.py. Collected across samples, so the report waits for them all.
    ch_report_files = ch_measurements
        .mix(
            QC_GATES.out.qc.map { meta, f -> [meta, 'qc', f] },
            READ_QC.out.gc_hist.map { meta, f -> [meta, 'gc-hist', f] },
            READ_QC.out.genomescope.flatMap { meta, fs ->
                (fs instanceof List ? fs : [fs]).findAll { it.name.endsWith('_linear_plot.png') && !it.name.contains('transformed') }.collect { f -> [meta, 'image-genomescope', f] }
            },
            ASSEMBLY.out.attempts.map { meta, f -> [meta, 'assembly-attempts', f] },
            ASSEMBLY.out.cluster_dotplots.flatMap { meta, fs ->
                (fs instanceof List ? fs : [fs]).collect { f -> [meta, 'image-dotplot', f] }
            },
            ASSEMBLY.out.autocycler_table.map { meta, f -> [meta, 'autocycler', f] },
            SCORING.out.scores.map { meta, f -> [meta, 'full-assemblies', f] },
            CHECKS.out.report,
        )
        .map { meta, kind, f -> [[meta.id, kind], f] }
        .toList()
        .filter { it }
        .multiMap { entries ->
            entries: entries.collect { it[0] }
            files: entries.collect { it[1] }
        }

    // Collated once every other process has run, and kept with the provenance reports.
    ch_versions = channel.topic('versions')
        .map { process, tool, version -> "${process.tokenize(':').last()}\t${tool}\t${version}" }
        .unique()
        .collectFile(name: 'software_versions.tsv', sort: true, newLine: true, storeDir: "${params.outdir}/pipeline_info")
    ch_run_info = channel.of(runInfoJson()).collectFile(name: 'run_info.json')

    COLLECT_METRICS(ch_report_files.entries, ch_report_files.files, ch_versions, ch_run_info)
    QUARTO_REPORT(COLLECT_METRICS.out.summary, files("${projectDir}/assets/report/*"), params.per_sample_reports)

    emit:
    reads = REMOVE_HUMAN.out.reads
    assembly = FINISHING.out.assembly
    contigs = FINISHING.out.contigs
    plasmid_audit = FINISHING.out.plasmid_audit
    assembly_source = SCORING.out.assembly_source
    read_qc = READ_QC.out.summary
    contamination = CONTAMINATION.out.summary
    qc = QC_GATES.out.qc
    report = QUARTO_REPORT.out.html
}

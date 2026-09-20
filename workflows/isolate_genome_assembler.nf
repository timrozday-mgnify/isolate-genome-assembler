include { READ_QC       } from '../subworkflows/local/read_qc'
include { CONTAMINATION } from '../subworkflows/local/contamination'
include { REMOVE_HUMAN  } from '../modules/local/remove_human'

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

    emit:
    reads = REMOVE_HUMAN.out.reads
    read_qc = READ_QC.out.summary
    contamination = CONTAMINATION.out.summary
}

// Collapse every read QC measurement for a sample into one tidy row, and record which
// genome size estimate was used. No thresholds are applied here: gating is `bin/qc_gates.py`.
process READ_QC_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(normalisation), path(seqkit_stats), path(gc_hist), path(duplicates), path(adapters), path(genome_size_autocycler), path(genomescope_summary)

    output:
    tuple val(meta), path("${meta.id}.read_qc.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def declared_genome_size = meta.genome_size ? "--declared-genome-size ${meta.genome_size}" : ''
    """
    read_qc_summary.py \\
        --sample ${meta.id} \\
        --normalisation ${normalisation} \\
        --seqkit-stats ${seqkit_stats} \\
        --gc-hist ${gc_hist} \\
        --duplicates ${duplicates} \\
        --adapters ${adapters} \\
        --genome-size-autocycler ${genome_size_autocycler} \\
        --genomescope-summary ${genomescope_summary} \\
        ${declared_genome_size} \\
        --output ${meta.id}.read_qc.tsv
    """

    stub:
    def declared_genome_size = meta.genome_size ? "--declared-genome-size ${meta.genome_size}" : ''
    """
    read_qc_summary.py \\
        --sample ${meta.id} \\
        --genome-size-autocycler ${genome_size_autocycler} \\
        ${declared_genome_size} \\
        --output ${meta.id}.read_qc.tsv
    """
}

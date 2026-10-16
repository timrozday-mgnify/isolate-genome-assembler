// Gene-level integrity (check E): the IDEEL protein length ratio and depth over each rRNA
// locus, from Bakta's annotation, the DIAMOND hits and mosdepth's windows.
process GENE_CHECKS {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(proteins), path(hits), path(gff), path(windows)  // windows optional

    output:
    tuple val(meta), path("${meta.id}.ideel.tsv"), emit: ratios
    tuple val(meta), path("${meta.id}.ideel_summary.tsv"), emit: summary
    tuple val(meta), path("${meta.id}.rrna_depth.tsv"), emit: rrna
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Without mapped reads there are no coverage windows, so rrna_depth.tsv comes out
    // empty and the IDEEL ratios, which need none, are still measured.
    def windows_arg = windows ? "--windows ${windows}" : ''
    """
    ideel_ratios.py \\
        --sample ${meta.id} \\
        --proteins ${proteins} \\
        --hits ${hits} \\
        --gff ${gff} \\
        ${windows_arg} \\
        --min-ratio ${params.ideel_min_ratio} \\
        --output ${meta.id}.ideel.tsv \\
        --summary ${meta.id}.ideel_summary.tsv \\
        --rrna-output ${meta.id}.rrna_depth.tsv
    """

    stub:
    def windows_stub = windows ? "--windows ${windows}" : ''
    """
    ideel_ratios.py \\
        --sample ${meta.id} \\
        --proteins ${proteins} \\
        --hits ${hits} \\
        --gff ${gff} \\
        ${windows_stub} \\
        --output ${meta.id}.ideel.tsv \\
        --summary ${meta.id}.ideel_summary.tsv \\
        --rrna-output ${meta.id}.rrna_depth.tsv
    """
}

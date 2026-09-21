// Bin read-vs-assembly disagreements by allele frequency, with homopolymer indels tallied
// separately (check D).
process VARIANT_SCAN {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(vcf), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.variants.tsv"), emit: variants
    tuple val(meta), path("${meta.id}.variant_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    variant_scan.py \\
        --sample ${meta.id} \\
        --vcf ${vcf} \\
        --assembly ${assembly} \\
        --min-depth ${params.variant_min_depth} \\
        --homopolymer ${params.homopolymer_min_len} \\
        --output ${meta.id}.variants.tsv \\
        --summary ${meta.id}.variant_summary.tsv
    """

    stub:
    """
    variant_scan.py \\
        --sample ${meta.id} \\
        --vcf ${vcf} \\
        --assembly ${assembly} \\
        --output ${meta.id}.variants.tsv \\
        --summary ${meta.id}.variant_summary.tsv
    """
}

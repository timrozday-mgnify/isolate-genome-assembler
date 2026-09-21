// Turn mosdepth's fixed windows into flagged low/high-depth regions and per-replicon depth
// relative to the chromosome (checks B and G).
process COVERAGE_REGIONS {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(windows), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.coverage_regions.tsv"), emit: regions
    tuple val(meta), path("${meta.id}.contig_depth.tsv"), emit: depth
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    coverage_regions.py \\
        --sample ${meta.id} \\
        --windows ${windows} \\
        --contigs ${contigs} \\
        --low ${params.coverage_low_ratio} \\
        --high ${params.coverage_high_ratio} \\
        --regions-output ${meta.id}.coverage_regions.tsv \\
        --depth-output ${meta.id}.contig_depth.tsv
    """

    stub:
    """
    coverage_regions.py \\
        --sample ${meta.id} \\
        --windows ${windows} \\
        --contigs ${contigs} \\
        --regions-output ${meta.id}.coverage_regions.tsv \\
        --depth-output ${meta.id}.contig_depth.tsv
    """
}

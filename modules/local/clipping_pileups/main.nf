// Positions where many reads are clipped at once: candidate misjoins (check B).
process CLIPPING_PILEUPS {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0'
        : 'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.clipping.tsv"), emit: pileups
    tuple val("${task.process}"), val('pysam'), eval('python3 -c "import pysam; print(pysam.__version__)"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    clipping_pileups.py \\
        --sample ${meta.id} \\
        --bam ${bam} \\
        --min-clip ${params.clip_min_length} \\
        --min-reads ${params.clip_min_reads} \\
        --min-fraction ${params.clip_min_fraction} \\
        --output ${meta.id}.clipping.tsv
    """

    stub:
    """
    printf 'sample\\tcontig\\tstart\\tend\\tclipped_reads\\tdepth\\tclipped_fraction\\n' > ${meta.id}.clipping.tsv
    """
}

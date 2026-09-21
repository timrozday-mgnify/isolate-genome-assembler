// Name the finished replicons and record what is wrong with each of them. Nothing is
// dropped unless --drop_flagged_contigs says so, and then it goes to removed_contigs.fasta.
process CLASSIFY_REPLICONS {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(rotated), path(unrotated), path(end_overlaps)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: assembly
    tuple val(meta), path("${meta.id}.contigs.tsv"), emit: table
    tuple val(meta), path("${meta.id}.removed_contigs.fasta"), emit: removed
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def unrotated_arg = unrotated ? "--unrotated ${unrotated}" : ''
    def drop = params.drop_flagged_contigs ? '--drop-flagged' : ''
    """
    classify_replicons.py \\
        --sample ${meta.id} \\
        --rotated ${rotated} \\
        ${unrotated_arg} \\
        --end-overlaps ${end_overlaps} \\
        --min-contig-len ${params.min_contig_len} \\
        --chromosome-min-len ${params.chromosome_min_len} \\
        --min-depth-ratio ${params.min_contig_depth_ratio} \\
        ${drop} \\
        --output ${meta.id}.fasta \\
        --removed ${meta.id}.removed_contigs.fasta \\
        --table ${meta.id}.contigs.tsv
    """

    stub:
    def unrotated_arg = unrotated ? "--unrotated ${unrotated}" : ''
    """
    classify_replicons.py \\
        --sample ${meta.id} \\
        --rotated ${rotated} \\
        ${unrotated_arg} \\
        --end-overlaps ${end_overlaps} \\
        --output ${meta.id}.fasta \\
        --removed ${meta.id}.removed_contigs.fasta \\
        --table ${meta.id}.contigs.tsv
    """
}

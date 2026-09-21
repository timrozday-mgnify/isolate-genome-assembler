// Which assemblers' contigs reached each QC-pass Autocycler cluster, and which of them
// `trim` then excluded (check C).
process ASSEMBLER_CONTRIBUTION {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(autocycler_dir, stageAs: 'autocycler_out')

    output:
    tuple val(meta), path("${meta.id}.assembler_contribution.tsv"), emit: contribution
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    assembler_contribution.py \\
        --sample ${meta.id} \\
        --autocycler-dir ${autocycler_dir} \\
        --output ${meta.id}.assembler_contribution.tsv
    """

    stub:
    """
    assembler_contribution.py \\
        --sample ${meta.id} \\
        --autocycler-dir ${autocycler_dir} \\
        --output ${meta.id}.assembler_contribution.tsv
    """
}

// Cut the duplicated sequence a circular contig carries at its ends, the one job of
// `autocycler trim` that a single FASTA can have done to it. Circularity comes out of the
// sequence here, not out of the assembler's tag, so every candidate is judged the same way.
process CIRCULARISE {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(assembly), path(overlaps)

    output:
    tuple val(meta), path("${meta.id}.circularised.fasta"), emit: assembly
    tuple val(meta), path("${meta.id}.circularity.tsv"), emit: circularity
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    circularise.py \\
        --assembly ${assembly} \\
        --overlaps ${overlaps} \\
        --min-identity ${params.circularise_min_identity} \\
        --output ${meta.id}.circularised.fasta \\
        --table ${meta.id}.circularity.tsv
    """

    stub:
    """
    cp ${assembly} ${meta.id}.circularised.fasta
    printf 'contig\\tlength_before\\toverlap\\tlength_after\\tidentity\\tcircular\\tflag\\n' > ${meta.id}.circularity.tsv
    """
}

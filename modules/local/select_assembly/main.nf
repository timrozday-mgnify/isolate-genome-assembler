// Pick the Autocycler consensus when it is fully resolved, and the full-read Flye assembly
// otherwise. The choice and its reason go into the sample's metrics; the warn status for a
// fallback is applied later, by bin/qc_gates.py.
process SELECT_ASSEMBLY {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(consensus), path(metrics), path(fallback)

    output:
    tuple val(meta), path("${meta.id}.assembly.fasta"), emit: assembly
    tuple val(meta), path("${meta.id}.assembly_source.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fallback_arg = fallback ? "--fallback ${fallback}" : ''
    """
    select_assembly.py \\
        --sample ${meta.id} \\
        --consensus ${consensus} \\
        --metrics ${metrics} \\
        ${fallback_arg} \\
        --output ${meta.id}.assembly.fasta \\
        --summary ${meta.id}.assembly_source.tsv
    """

    stub:
    def fallback_arg = fallback ? "--fallback ${fallback}" : ''
    """
    select_assembly.py \\
        --sample ${meta.id} \\
        --consensus ${consensus} \\
        --metrics ${metrics} \\
        ${fallback_arg} \\
        --output ${meta.id}.assembly.fasta \\
        --summary ${meta.id}.assembly_source.tsv
    """
}

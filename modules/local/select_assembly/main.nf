// Pick which candidate assembly the sample is finished from: the Autocycler consensus, or
// the best-scoring full-read assembly, as --assembly_selection says. The choice and its
// reason go into the sample's metrics; the warn status for a fallback is applied later, by
// bin/qc_gates.py.
process SELECT_ASSEMBLY {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(consensus), path(metrics), path(fallback), path(scores), val(assemblers), path(candidates, stageAs: 'candidates/?/*')

    output:
    tuple val(meta), path("${meta.id}.assembly.fasta"), emit: assembly
    tuple val(meta), path("${meta.id}.assembly_source.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fallback_arg = fallback ? "--fallback ${fallback}" : ''
    def scores_arg = scores ? "--scores ${scores}" : ''
    def files = candidates instanceof List ? candidates : [candidates]
    def named = assemblers && candidates
        ? "--candidates " + [assemblers, files].transpose().collect { assembler, file -> "${assembler}=${file}" }.join(' ')
        : ''
    """
    select_assembly.py \\
        --sample ${meta.id} \\
        --consensus ${consensus} \\
        --metrics ${metrics} \\
        --selection ${params.assembly_selection} \\
        ${fallback_arg} \\
        ${scores_arg} \\
        ${named} \\
        --output ${meta.id}.assembly.fasta \\
        --summary ${meta.id}.assembly_source.tsv
    """

    stub:
    def fallback_arg = fallback ? "--fallback ${fallback}" : ''
    def scores_arg = scores ? "--scores ${scores}" : ''
    def files = candidates instanceof List ? candidates : [candidates]
    def named = assemblers && candidates
        ? "--candidates " + [assemblers, files].transpose().collect { assembler, file -> "${assembler}=${file}" }.join(' ')
        : ''
    """
    select_assembly.py \\
        --sample ${meta.id} \\
        --consensus ${consensus} \\
        --metrics ${metrics} \\
        --selection ${params.assembly_selection} \\
        ${fallback_arg} \\
        ${scores_arg} \\
        ${named} \\
        --output ${meta.id}.assembly.fasta \\
        --summary ${meta.id}.assembly_source.tsv
    """
}

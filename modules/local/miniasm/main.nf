// miniasm's unpolished overlap-layout graph. It has no consensus step of its own, so
// MINIPOLISH follows.
process MINIASM {
    tag "${meta.id}_${subset}"
    label 'process_medium'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/miniasm:0.3--h577a1d6_5'
        : 'quay.io/biocontainers/miniasm:0.3--h577a1d6_5'}"

    input:
    tuple val(meta), val(subset), path(reads), path(overlap)

    output:
    tuple val(meta), val(subset), path(reads), path("${meta.id}_${subset}.unpolished.gfa"), emit: graph
    tuple val("${task.process}"), val('miniasm'), eval("miniasm -V"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    miniasm -f ${reads} ${args} ${overlap} > ${meta.id}_${subset}.unpolished.gfa
    """

    stub:
    """
    touch ${meta.id}_${subset}.unpolished.gfa
    """
}

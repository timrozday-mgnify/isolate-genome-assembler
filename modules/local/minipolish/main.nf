// Racon polishing of the miniasm graph, which also puts read depths on the segments.
process MINIPOLISH {
    tag "${meta.id}_${subset}"
    label 'process_medium'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/minipolish:0.2.1--pyhdfd78af_0'
        : 'quay.io/biocontainers/minipolish:0.2.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), val(subset), path(reads), path(graph)

    output:
    tuple val(meta), val('miniasm'), val(subset), path('out/*'), emit: assembly
    tuple val("${task.process}"), val('minipolish'), eval("minipolish --version | sed 's/^.*v//'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p out
    minipolish \\
        --threads ${task.cpus} \\
        --minimap2-preset map-hifi \\
        ${args} \\
        ${reads} \\
        ${graph} \\
        > out/miniasm.gfa
    """

    stub:
    """
    mkdir -p out
    printf 'S\\tutg1c\\tACGT\\tdp:f:30\\n' > out/miniasm.gfa
    """
}

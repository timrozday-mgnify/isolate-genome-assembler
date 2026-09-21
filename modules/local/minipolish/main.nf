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
    // --skip_initial: the per-segment round aborts the whole run when Racon returns nothing
    // for a single-read segment (Minipolish 0.2.1 exits before its own drop-the-segment
    // fallback). The full rounds still polish every segment against all the reads.
    """
    mkdir -p out
    minipolish \\
        --threads ${task.cpus} \\
        --minimap2-preset map-hifi \\
        --skip_initial \\
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

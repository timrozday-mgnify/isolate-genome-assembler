// All-vs-all read overlap for the miniasm arm. The HiFi preset is Autocycler helper.rs's
// `-k23 -Xw11 -e0 -m100` rather than the ava-pb preset, which is tuned for noisy CLR reads.
process MINIASM_OVERLAP {
    tag "${meta.id}_${subset}"
    label 'process_medium'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/minipolish:0.2.1--pyhdfd78af_0'
        : 'quay.io/biocontainers/minipolish:0.2.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), val(subset), path(reads)

    output:
    tuple val(meta), val(subset), path(reads), path("${meta.id}_${subset}.overlap.paf"), emit: overlap
    tuple val("${task.process}"), val('minimap2'), eval("minimap2 --version"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-k23 -Xw11 -e0 -m100'
    """
    minimap2 -t ${task.cpus} ${args} ${reads} ${reads} > ${meta.id}_${subset}.overlap.paf
    """

    stub:
    """
    touch ${meta.id}_${subset}.overlap.paf
    """
}

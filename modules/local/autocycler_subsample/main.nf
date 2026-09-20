// Split the read set into independent subsets, each deep enough to assemble on its own.
// Independent subsets are what makes the consensus meaningful: errors that are random
// per-assembly get outvoted.
process AUTOCYCLER_SUBSAMPLE {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/autocycler:0.7.0--h79ce301_0'
        : 'quay.io/biocontainers/autocycler:0.7.0--h79ce301_0'}"

    input:
    tuple val(meta), path(reads), val(genome_size)

    output:
    tuple val(meta), val(genome_size), path("subsets/*.fastq"), emit: subsets
    tuple val(meta), path("${meta.id}.subsample.yaml"), emit: metrics
    tuple val("${task.process}"), val('autocycler'), eval("autocycler --version | sed 's/^.*r //; s/^autocycler //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    autocycler subsample \\
        --reads ${reads} \\
        --out_dir subsets \\
        --genome_size ${genome_size} \\
        --count ${params.subsample_count} \\
        --min_read_depth ${params.subsample_min_depth} \\
        ${args}

    cp subsets/subsample.yaml ${meta.id}.subsample.yaml
    """

    stub:
    """
    mkdir -p subsets
    for subset in \$(seq 1 ${params.subsample_count}); do
        touch "subsets/sample_\$(printf '%02d' "\$subset").fastq"
    done
    echo 'input_read_count: 1000' > ${meta.id}.subsample.yaml
    """
}

// Autocycler's own genome size estimate (a Raven assembly of all reads). Used as the
// primary estimate when the samplesheet does not give `genome_size`, and cross-checked
// against the k-mer estimate from KMC + GenomeScope2.
process AUTOCYCLER_GENOME_SIZE {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/autocycler:0.7.0--h79ce301_0'
        : 'quay.io/biocontainers/autocycler:0.7.0--h79ce301_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.genome_size_autocycler.txt"), emit: genome_size
    tuple val("${task.process}"), val('autocycler'), eval("autocycler --version | sed 's/^.*r //; s/^autocycler //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    autocycler helper genome_size \\
        --reads ${reads} \\
        --threads ${task.cpus} \\
        ${args} \\
        > ${meta.id}.genome_size_autocycler.txt
    """

    stub:
    """
    echo 5000000 > ${meta.id}.genome_size_autocycler.txt
    """
}

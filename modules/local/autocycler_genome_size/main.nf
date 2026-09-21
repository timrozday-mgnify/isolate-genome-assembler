// Autocycler's own genome size estimate (a Raven assembly of all reads). Used as the
// primary estimate when the samplesheet does not give `genome_size`, and cross-checked
// against the k-mer estimate from KMC + GenomeScope2.
process AUTOCYCLER_GENOME_SIZE {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/raven-assembler:1.8.3--h5ca1c30_3'
        : 'quay.io/biocontainers/raven-assembler:1.8.3--h5ca1c30_3'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.genome_size_autocycler.txt"), emit: genome_size
    tuple val("${task.process}"), val('raven'), eval("raven --version"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // ponytail: same as `autocycler helper genome_size` (helper.rs genome_size_raven), which needs raven on PATH
    // and the autocycler biocontainer does not have it.
    """
    raven --threads ${task.cpus} --disable-checkpoints --polishing-rounds 1 ${args} ${reads} > assembly.fasta
    awk '!/^>/ { total += length(\$0) } END { if (total == 0) exit 1; print total }' assembly.fasta \\
        > ${meta.id}.genome_size_autocycler.txt
    rm assembly.fasta
    """

    stub:
    """
    echo 5000000 > ${meta.id}.genome_size_autocycler.txt
    """
}

// Raven's overlap-layout assembly. As in Autocycler's helper.rs the contigs are taken from
// stdout, with the GFA kept alongside for the report.
process RAVEN {
    tag "${meta.id}_${subset}"
    label 'process_medium'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/raven-assembler:1.8.3--h5ca1c30_3'
        : 'quay.io/biocontainers/raven-assembler:1.8.3--h5ca1c30_3'}"

    input:
    tuple val(meta), val(subset), path(reads)

    output:
    tuple val(meta), val('raven'), val(subset), path('out/*'), emit: assembly
    tuple val("${task.process}"), val('raven'), eval("raven --version"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p out
    raven \\
        --threads ${task.cpus} \\
        --disable-checkpoints \\
        --graphical-fragment-assembly out/raven.gfa \\
        ${args} \\
        ${reads} \\
        > out/raven.fasta
    """

    stub:
    """
    mkdir -p out
    printf '>Utg1\\nACGT\\n' > out/raven.fasta
    printf 'S\\tUtg1\\tACGT\\n' > out/raven.gfa
    """
}

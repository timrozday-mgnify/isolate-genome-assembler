// LJA's multiplex de Bruijn graph, built for HiFi reads. Off by default; flags follow
// Autocycler v0.7.0's helper.rs.
process LJA {
    tag "${meta.id}_${subset}"
    label 'process_high'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/lja:0.2--h5b5514e_2'
        : 'quay.io/biocontainers/lja:0.2--h5b5514e_2'}"

    input:
    tuple val(meta), val(subset), path(reads)

    output:
    tuple val(meta), val('lja'), val(subset), path('out/*'), emit: assembly
    // lja has no --version flag; the version is the container's.
    tuple val("${task.process}"), val('lja'), val('0.2'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    lja --output-dir lja_out --reads ${reads} --threads ${task.cpus} ${args}
    mkdir -p out
    cp lja_out/assembly.fasta out/
    """

    stub:
    """
    mkdir -p out
    printf '>1\\nACGT\\n' > out/assembly.fasta
    """
}

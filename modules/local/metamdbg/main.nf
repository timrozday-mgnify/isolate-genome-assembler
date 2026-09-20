// metaMDBG's minimizer-space de Bruijn graph: the algorithmically most distinct member of
// the set, which is what makes its errors unlikely to coincide with the others'.
process METAMDBG {
    tag "${meta.id}_${subset}"
    label 'process_medium'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/metamdbg:1.4--h3be2455_0'
        : 'quay.io/biocontainers/metamdbg:1.4--h3be2455_0'}"

    input:
    tuple val(meta), val(subset), path(reads)

    output:
    tuple val(meta), val('metamdbg'), val(subset), path('out/*'), emit: assembly
    tuple val("${task.process}"), val('metaMDBG'), eval("metaMDBG --version 2>&1 | head -n 1"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    metaMDBG asm --out-dir out --in-hifi ${reads} --threads ${task.cpus} ${args}
    find out -mindepth 1 -maxdepth 1 ! -name 'contigs.fasta.gz' ! -name 'metaMDBG.log' \\
        -exec rm -rf {} +
    """

    stub:
    """
    mkdir -p out
    printf '>ctg1\\nACGT\\n' | gzip -c > out/contigs.fasta.gz
    """
}

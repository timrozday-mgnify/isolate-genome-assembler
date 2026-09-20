// K-mer count histogram, the input GenomeScope2 fits for an independent genome size
// estimate and for the heterozygosity/mixed-strain peak.
process KMC_HISTOGRAM {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/kmc:3.2.4--h5ca1c30_4'
        : 'quay.io/biocontainers/kmc:3.2.4--h5ca1c30_4'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.kmc.histo"), emit: histogram
    tuple val("${task.process}"), val('kmc'), eval("kmc 2>&1 | sed -n 's/^K-Mer Counter (KMC) ver. \\\\([0-9.]*\\\\).*/\\\\1/p'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "-k${params.kmer_size} -ci1 -cs10000"
    """
    mkdir -p kmc_tmp
    kmc ${args} -t${task.cpus} -m${(task.memory.toGiga() as int)} -fq ${reads} ${meta.id}.kmc kmc_tmp
    kmc_tools transform ${meta.id}.kmc histogram ${meta.id}.kmc.histo -cx10000
    """

    stub:
    """
    printf '1\\t1000\\n2\\t500\\n' > ${meta.id}.kmc.histo
    """
}

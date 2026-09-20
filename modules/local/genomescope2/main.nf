// GenomeScope2 fit of the k-mer histogram. `-p 1` because an isolate is haploid; a second
// peak at half depth then means mixed strains rather than heterozygosity.
process GENOMESCOPE2 {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/genomescope2:2.1.0--py313r44hdfd78af_0'
        : 'quay.io/biocontainers/genomescope2:2.1.0--py313r44hdfd78af_0'}"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("genomescope/${meta.id}_summary.txt"), emit: summary
    tuple val(meta), path("genomescope/*"), emit: results
    tuple val("${task.process}"), val('genomescope2'), val('2.1.0'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-p 1'
    """
    genomescope2 \\
        --input ${histogram} \\
        --kmer_length ${params.kmer_size} \\
        --output genomescope \\
        --name_prefix ${meta.id} \\
        ${args}
    """

    stub:
    """
    mkdir -p genomescope
    printf 'GenomeScope version 2.0\\nGenome Haploid Length 5,000,000 bp\\n' > genomescope/${meta.id}_summary.txt
    """
}

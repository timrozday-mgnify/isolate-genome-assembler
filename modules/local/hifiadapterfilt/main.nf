// Screen HiFi reads for leftover PacBio adapter/primer sequence. The filtered read set is
// used downstream only when --remove_adapter_reads is true; the counts are always reported.
process HIFIADAPTERFILT {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/hifiadapterfilt:3.0.0--hdfd78af_0'
        : 'quay.io/biocontainers/hifiadapterfilt:3.0.0--hdfd78af_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filt.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.adapters.tsv"), emit: stats
    tuple val("${task.process}"), val('HiFiAdapterFilt'), val('3.0.0'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ln -s ${reads} ${meta.id}.fastq.gz
    pbadapterfilt.sh -p ${meta.id} -t ${task.cpus} ${args}

    # The report file records total reads, reads with adapter and the percentage.
    awk 'BEGIN { print "reads\\tadapter_reads\\tadapter_fraction" }
         /Number of ccs reads:/ { total = \$NF }
         /Number of adapter contaminated ccs reads:/ { contaminated = \$(NF - 2) }
         END { print total "\\t" contaminated "\\t" (total > 0 ? contaminated / total : 0) }' \\
        ${meta.id}.stats > ${meta.id}.adapters.tsv
    """

    stub:
    """
    echo | gzip -c > ${meta.id}.filt.fastq.gz
    printf 'reads\\tadapter_reads\\tadapter_fraction\\n1000\\t0\\t0\\n' > ${meta.id}.adapters.tsv
    """
}

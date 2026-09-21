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
    // Thresholds for the NGB00972 adapter; NGB00973 is fixed at 97% over 34 bp, as in pbadapterfilt.sh.
    def min_length = task.ext.min_length ?: 44
    def min_identity = task.ext.min_identity ?: 97
    // ponytail: runs pbadapterfilt.sh's steps directly. The 3.0.0 script misreads the work dir, relies on
    // GNU sed/dirname that the BusyBox container lacks, and cannot find its own BLAST DB here.
    """
    db=\$(dirname \$(command -v pbadapterfilt.sh))/DB/pacbio_vectors_db

    zcat -f ${reads} | awk 'NR % 4 == 1 { print ">" substr(\$1, 2) } NR % 4 == 2' > reads.fasta
    blastn -db \$db -query reads.fasta -num_threads ${task.cpus} -task blastn -reward 1 -penalty -5 \\
        -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue 700 -searchsp 1750000000000 \\
        -outfmt 6 > ${meta.id}.blastout
    awk -v len=${min_length} -v pct=${min_identity} \\
        '(\$2 ~ /NGB00972/ && \$3 >= pct && \$4 >= len) || (\$2 ~ /NGB00973/ && \$3 >= 97 && \$4 >= 34) { print \$1 }' \\
        ${meta.id}.blastout | sort -u > ${meta.id}.blocklist

    zcat -f ${reads} \\
        | awk 'FILENAME == "${meta.id}.blocklist" { block[\$1]; next }
               FNR % 4 == 1 { keep = !(substr(\$1, 2) in block) }
               keep' ${meta.id}.blocklist - \\
        | gzip -1 > ${meta.id}.filt.fastq.gz

    reads=\$(grep -c '^>' reads.fasta || true)
    contaminated=\$(wc -l < ${meta.id}.blocklist)
    printf 'reads\\tadapter_reads\\tadapter_fraction\\n' > ${meta.id}.adapters.tsv
    awk -v n=\$reads -v c=\$contaminated 'BEGIN { print n "\\t" c "\\t" (n > 0 ? c / n : 0) }' >> ${meta.id}.adapters.tsv
    rm reads.fasta
    """

    stub:
    """
    echo | gzip -c > ${meta.id}.filt.fastq.gz
    printf 'reads\\tadapter_reads\\tadapter_fraction\\n1000\\t0\\t0\\n' > ${meta.id}.adapters.tsv
    """
}

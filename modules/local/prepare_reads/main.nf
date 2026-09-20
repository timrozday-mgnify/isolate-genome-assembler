// Normalise a sample's HiFi inputs into one FASTQ: BAM -> FASTQ, concatenate runs, and
// drop reads whose mean accuracy is below --min_read_qv. A hifi_reads file should not
// contain any, so the count of dropped reads is itself a QC signal.
process PREPARE_READS {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_1'
        : 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.reads.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.read_normalisation.tsv"), emit: stats
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # ponytail: one awk pass does the QV filter for both sources; a Python script here
    # would need a Python image just to average quality strings.
    for input in ${reads}; do
        case "\$input" in
            *.bam) samtools fastq -@ ${task.cpus - 1} "\$input" ;;
            *.gz)  gzip -cd "\$input" ;;
            *)     cat "\$input" ;;
        esac
    done \\
        | awk -v min_qv=${params.min_read_qv} -v stats='${meta.id}.read_normalisation.tsv' '
            BEGIN { for (i = 33; i < 127; i++) phred[sprintf("%c", i)] = i - 33 }
            NR % 4 == 1 { header = \$0 }
            NR % 4 == 2 { sequence = \$0 }
            NR % 4 == 0 {
                total++
                length_qual = length(\$0)
                error_sum = 0
                for (i = 1; i <= length_qual; i++) {
                    error_sum += 10 ^ (-phred[substr(\$0, i, 1)] / 10)
                }
                read_qv = length_qual > 0 ? -10 * log(error_sum / length_qual) / log(10) : 0
                if (read_qv < min_qv) { low_qv++; next }
                kept++
                bases += length(sequence)
                print header "\\n" sequence "\\n+\\n" \$0
            }
            END {
                printf "reads_in\\treads_kept\\treads_below_min_qv\\tbases_kept\\tmin_read_qv\\n" > stats
                printf "%d\\t%d\\t%d\\t%d\\t%s\\n", total, kept, low_qv, bases, min_qv > stats
            }' \\
        | gzip -c > ${meta.id}.reads.fastq.gz
    """

    stub:
    """
    echo | gzip -c > ${meta.id}.reads.fastq.gz
    printf 'reads_in\\treads_kept\\treads_below_min_qv\\tbases_kept\\tmin_read_qv\\n1000\\t1000\\t0\\t150000000\\t${params.min_read_qv}\\n' > ${meta.id}.read_normalisation.tsv
    """
}

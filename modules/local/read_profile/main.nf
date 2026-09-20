// Database-free read profile: per-read GC histogram (a wide or bimodal distribution is an
// early contamination signal) and the exact-duplicate fraction by id and by sequence.
process READ_PROFILE {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--he881be0_0'
        : 'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.gc_hist.tsv"), emit: gc_hist
    tuple val(meta), path("${meta.id}.duplicates.tsv"), emit: duplicates
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    seqkit fx2tab --name --only-id --gc --threads ${task.cpus} ${reads} \\
        | awk 'BEGIN { print "gc_percent_bin\\tread_count" }
               { bin[int(\$2)]++ }
               END { for (gc = 0; gc <= 100; gc++) if (bin[gc]) print gc "\\t" bin[gc] }' \\
        > ${meta.id}.gc_hist.tsv

    total=\$(seqkit stats --tabular ${reads} | awk 'NR == 2 { print \$4 }')
    unique_ids=\$(seqkit rmdup --by-name --threads ${task.cpus} ${reads} 2> /dev/null | seqkit stats --tabular | awk 'NR == 2 { print \$4 }')
    unique_seqs=\$(seqkit rmdup --by-seq --threads ${task.cpus} ${reads} 2> /dev/null | seqkit stats --tabular | awk 'NR == 2 { print \$4 }')
    printf 'reads\\tduplicate_ids\\tduplicate_sequences\\n%s\\t%s\\t%s\\n' \\
        "\$total" "\$((total - unique_ids))" "\$((total - unique_seqs))" \\
        > ${meta.id}.duplicates.tsv
    """

    stub:
    """
    printf 'gc_percent_bin\\tread_count\\n50\\t1000\\n' > ${meta.id}.gc_hist.tsv
    printf 'reads\\tduplicate_ids\\tduplicate_sequences\\n1000\\t0\\t0\\n' > ${meta.id}.duplicates.tsv
    """
}

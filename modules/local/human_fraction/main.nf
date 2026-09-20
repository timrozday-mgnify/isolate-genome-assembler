// Fraction of reads with any alignment to the human reference, and their names so they can
// be removed. This is a direct measure, independent of sylph's containment ANI threshold.
process HUMAN_FRACTION {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--he881be0_0'
        : 'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'}"

    input:
    tuple val(meta), path(paf), path(reads)

    output:
    tuple val(meta), path("${meta.id}.human_fraction.tsv"), emit: fraction
    tuple val(meta), path("${meta.id}.human_read_ids.txt"), emit: read_ids
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cut -f1 ${paf} | sort --unique > ${meta.id}.human_read_ids.txt
    total=\$(seqkit stats --tabular ${reads} | awk 'NR == 2 { print \$4 }')
    human=\$(wc -l < ${meta.id}.human_read_ids.txt)
    printf 'reads\\thuman_reads\\thuman_fraction\\n%s\\t%s\\t%s\\n' \\
        "\$total" "\$human" "\$(awk -v h="\$human" -v t="\$total" 'BEGIN { print t > 0 ? h / t : 0 }')" \\
        > ${meta.id}.human_fraction.tsv
    """

    stub:
    """
    : > ${meta.id}.human_read_ids.txt
    printf 'reads\\thuman_reads\\thuman_fraction\\n1000\\t0\\t0\\n' > ${meta.id}.human_fraction.tsv
    """
}

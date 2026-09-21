// Self-align the ends of every contig. Autocycler's `trim` removes the duplicated sequence
// a circular contig carries at its ends; the Flye fallback path is the one that can still
// have it, and a contig that does is over-long by the size of that overlap.
process CONTIG_ENDS {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/minimap2:2.30--h577a1d6_0'
        : 'quay.io/biocontainers/minimap2:2.30--h577a1d6_0'}"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.end_overlaps.paf"), emit: overlaps
    tuple val("${task.process}"), val('minimap2'), eval("minimap2 --version"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-x asm5'
    """
    # The first and last ${params.contig_end_window} bases of each contig, as their own
    # sequences, so an end-to-end match shows up as a <contig>_start vs <contig>_end hit.
    awk -v window=${params.contig_end_window} '
        /^>/ {
            if (name != "") print_ends()
            name = substr(\$1, 2); seq = ""; next
        }
        { seq = seq \$0 }
        END { if (name != "") print_ends() }
        function print_ends(   len, head, tail) {
            len = length(seq)
            head = (len > window) ? substr(seq, 1, window) : seq
            tail = (len > window) ? substr(seq, len - window + 1) : seq
            print ">" name "_start\\n" head
            print ">" name "_end\\n" tail
        }' ${assembly} > ends.fasta

    minimap2 -t ${task.cpus} ${args} ends.fasta ends.fasta \\
        | awk -F'\\t' '\$1 != \$6' \\
        > ${meta.id}.end_overlaps.paf
    """

    stub:
    """
    touch ${meta.id}.end_overlaps.paf
    """
}

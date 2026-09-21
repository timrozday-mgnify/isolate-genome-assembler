// Inspector's read-based structural and small-scale error calls, and its QV (check B).
// Its defaults skip contigs under 10 kb, and structural calls on contigs under 1 Mb,
// which would exclude every plasmid; an isolate genome needs both looked at.
process INSPECTOR {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/inspector:1.3.1--hdfd78af_1'
        : 'quay.io/biocontainers/inspector:1.3.1--hdfd78af_1'}"

    input:
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path("${meta.id}.inspector"), emit: outdir
    tuple val(meta), path("${meta.id}.inspector/summary_statistics"), emit: summary
    // inspector.py --version still prints v1.0.1 in the 1.3.1 release, so the version is
    // pinned here with the container.
    tuple val("${task.process}"), val('inspector'), val('1.3.1'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    inspector.py \\
        -c ${assembly} \\
        -r ${reads} \\
        -d hifi \\
        -t ${task.cpus} \\
        -o ${meta.id}.inspector \\
        ${args}
    # The read alignments duplicate MAP_READS' BAM and are the bulk of the directory.
    rm -f ${meta.id}.inspector/read_to_contig.bam*
    """

    stub:
    """
    mkdir -p ${meta.id}.inspector
    printf 'Structural error\\t0\\n\\nQV\\t60\\n' > ${meta.id}.inspector/summary_statistics
    """
}

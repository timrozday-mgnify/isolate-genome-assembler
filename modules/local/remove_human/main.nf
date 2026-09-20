// Drop human reads when the contamination screen asked for it. The process always runs so
// the DAG stays static and resume-safe; the decision file says whether it does anything.
process REMOVE_HUMAN {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--he881be0_0'
        : 'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'}"

    input:
    tuple val(meta), path(reads), path(read_ids), path(decision)

    output:
    tuple val(meta), path("${meta.id}.clean.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [ "\$(cat ${decision})" = "true" ]; then
        seqkit grep --invert-match --by-name --pattern-file ${read_ids} --threads ${task.cpus} ${reads} \\
            | gzip -c > ${meta.id}.clean.fastq.gz
    else
        cp ${reads} ${meta.id}.clean.fastq.gz
    fi
    """

    stub:
    """
    echo | gzip -c > ${meta.id}.clean.fastq.gz
    """
}

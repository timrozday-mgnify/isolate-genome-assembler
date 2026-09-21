// Compare the plasmids Plassembler found in the full read set with the finished assembly.
// A plasmid that is missing is a warn, not a failure: the report shows the contig so a
// person can decide whether to add it back.
process PLASMID_AUDIT {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(plassembler_summary), path(skani)

    output:
    tuple val(meta), path("${meta.id}.plasmid_audit.tsv"), emit: audit
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def summary_arg = plassembler_summary ? "--plassembler-summary ${plassembler_summary}" : ''
    """
    plasmid_audit.py \\
        --sample ${meta.id} \\
        ${summary_arg} \\
        --skani ${skani} \\
        --min-identity ${params.plasmid_audit_min_identity} \\
        --min-coverage ${params.plasmid_audit_min_coverage} \\
        --output ${meta.id}.plasmid_audit.tsv
    """

    stub:
    def summary_arg = plassembler_summary ? "--plassembler-summary ${plassembler_summary}" : ''
    """
    plasmid_audit.py \\
        --sample ${meta.id} \\
        ${summary_arg} \\
        --skani ${skani} \\
        --output ${meta.id}.plasmid_audit.tsv
    """
}

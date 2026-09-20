// Combine the sylph GTDB profile, the human query and the human read fraction into one
// JSON of measurements, plus the `--remove_human auto` decision the read path acts on.
process CONTAMINATION_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(sylph_profile), path(sylph_tax), path(human_query), path(human_fraction)

    output:
    tuple val(meta), path("${meta.id}.contamination_summary.json"), emit: summary
    tuple val(meta), path("${meta.id}.remove_human.txt"), emit: remove_human
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def expected_taxon = meta.expected_taxon ? "--expected-taxon '${meta.expected_taxon}'" : ''
    """
    contamination_summary.py \\
        --sample ${meta.id} \\
        --sylph-profile ${sylph_profile} \\
        --sylph-tax ${sylph_tax} \\
        --human-query ${human_query} \\
        --human-fraction ${human_fraction} \\
        --remove-human ${params.remove_human} \\
        --max-human-fraction ${params.contam_max_human} \\
        ${expected_taxon} \\
        --output ${meta.id}.contamination_summary.json \\
        --remove-human-decision ${meta.id}.remove_human.txt
    """

    stub:
    """
    contamination_summary.py \\
        --sample ${meta.id} \\
        --remove-human ${params.remove_human} \\
        --max-human-fraction ${params.contam_max_human} \\
        --output ${meta.id}.contamination_summary.json \\
        --remove-human-decision ${meta.id}.remove_human.txt
    """
}

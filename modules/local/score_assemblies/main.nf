// Rank one sample's candidate assemblies -- every assembler's full-read assembly and the
// Autocycler consensus -- against each other, from checks that have already run.
process SCORE_ASSEMBLIES {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(stats, stageAs: 'stats/*'), path(circularity, stageAs: 'circularity/*'), path(qv, stageAs: 'qv/*'), path(completeness, stageAs: 'completeness/*'), path(mapping, stageAs: 'mapping/*'), path(clipping, stageAs: 'clipping/*'), val(genome_size)

    output:
    tuple val(meta), path("${meta.id}.full_assemblies.tsv"), emit: scores
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    score_assemblies.py \\
        --sample ${meta.id} \\
        --genome-size ${genome_size} \\
        --size-tolerance ${params.score_size_tolerance} \\
        --min-qv ${params.score_min_qv} \\
        --output ${meta.id}.full_assemblies.tsv
    """

    stub:
    """
    score_assemblies.py \\
        --sample ${meta.id} \\
        --genome-size ${genome_size} \\
        --output ${meta.id}.full_assemblies.tsv
    """
}

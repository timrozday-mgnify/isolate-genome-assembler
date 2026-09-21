// Stage 7: apply the thresholds to every measurement taken for one sample. The inputs
// arrive as parallel lists of names and files, because which checks produced a file
// differs between samples (a reference is optional; an unmapped-read assembly may fail).
process QC_GATES {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), val(names), path(metrics, stageAs: 'metrics/?/*')
    val thresholds // JSON text, converted from the thresholds YAML by the workflow

    output:
    tuple val(meta), path("${meta.id}.qc.json"), emit: qc
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def files = metrics instanceof List ? metrics : [metrics]
    def inputs = [names, files].transpose().collect { name, metric -> "--${name} ${metric}" }.join(' ')
    def expected = meta.expected_taxon ? "--expected-taxon '${meta.expected_taxon}'" : ''
    """
    cat > thresholds.json <<'JSON'
    ${thresholds}
    JSON

    qc_gates.py \\
        --sample ${meta.id} \\
        --thresholds thresholds.json \\
        ${expected} \\
        ${inputs} \\
        --output ${meta.id}.qc.json
    """

    stub:
    def files = metrics instanceof List ? metrics : [metrics]
    def inputs = [names, files].transpose().collect { name, metric -> "--${name} ${metric}" }.join(' ')
    """
    cat > thresholds.json <<'JSON'
    ${thresholds}
    JSON

    qc_gates.py \\
        --sample ${meta.id} \\
        --thresholds thresholds.json \\
        ${inputs} \\
        --output ${meta.id}.qc.json
    """
}

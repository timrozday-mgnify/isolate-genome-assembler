// Stage 8: every sample's measurements, verdicts and images as one run_summary/ bundle,
// the only thing the report reads. Files arrive as one list with a parallel list of
// [sample, kind] pairs, since which files exist differs between samples.
process COLLECT_METRICS {
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    val entries // [sample id, kind] per file
    path files, stageAs: 'inputs/?/*'
    path software_versions // process, tool, version TSV
    path run_info // JSON: pipeline, run and params

    output:
    // No versions topic output: this process reads the collated topic, so writing to it
    // would keep the topic open forever.
    path 'run_summary', emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def paths = files instanceof List ? files : [files]
    // Ids are validated in main.nf and kinds are fixed words, so neither needs quoting.
    def manifest = [entries, paths].transpose().collect { entry, f -> "${entry[0]} ${entry[1]} ${f}" }.join(' ')
    """
    printf '%s\\t%s\\t%s\\n' ${manifest} > manifest.tsv

    collect_metrics.py \\
        --manifest manifest.tsv \\
        --software-versions ${software_versions} \\
        --run-info ${run_info} \\
        --outdir run_summary
    """
}

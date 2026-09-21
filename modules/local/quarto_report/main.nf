// Stage 8: render the run report, plus one report per sample, from run_summary/ alone.
// embed-resources makes each a single self-contained HTML file.
process QUARTO_REPORT {
    label 'process_single'

    container 'ghcr.io/timrozday-mgnify/isolate-genome-assembler-report:latest'

    input:
    path run_summary
    path template, stageAs: 'template/*' // assets/report/*
    val per_sample

    // No versions topic output, for the same reason as COLLECT_METRICS: the report is
    // rendered after the topic is collated.
    output:
    path 'isolate_assembly_report.html', emit: html
    path 'samples/*.report.html', optional: true, emit: sample_html

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Quarto and Jupyter write caches under HOME, which is not writable in every container.
    export HOME=\$PWD XDG_CACHE_HOME=\$PWD/.cache XDG_DATA_HOME=\$PWD/.local/share

    # Quarto writes next to the .qmd, so render from a real copy, not the staged links.
    cp -rL template render
    summary=\$PWD/${run_summary}

    quarto render render/isolate_assembly_report.qmd -P summary_dir:\$summary
    mv render/isolate_assembly_report.html .

    mkdir samples
    if [ "${per_sample}" = true ]; then
        for id in \$(tail -n +2 ${run_summary}/samples.tsv | cut -f 1); do
            # --output breaks embed-resources (Quarto looks for the libs under the new
            # name), so each render takes the default name and is renamed.
            quarto render render/isolate_assembly_report.qmd -P summary_dir:\$summary -P sample:\$id
            mv render/isolate_assembly_report.html samples/\$id.report.html
        done
    fi
    """

    stub:
    """
    touch isolate_assembly_report.html
    mkdir samples
    if [ "${per_sample}" = true ]; then
        for id in \$(tail -n +2 ${run_summary}/samples.tsv | cut -f 1); do
            touch samples/\$id.report.html
        done
    fi
    """
}

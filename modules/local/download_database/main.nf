// Fetch one database file. `--prepare_databases` runs this once per site; the normal
// pipeline only reads the resulting paths through the database params.
process DOWNLOAD_DATABASE {
    tag "${name}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--hb829ee6_10'
        : 'quay.io/biocontainers/gnu-wget:1.18--hb829ee6_10'}"

    input:
    tuple val(name), val(url), val(filename)

    output:
    tuple val(name), path(filename), emit: database
    tuple val("${task.process}"), val('wget'), eval("wget --version | sed -n '1s/^GNU Wget \\\\([0-9.]*\\\\).*/\\\\1/p'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget --no-verbose --continue --output-document '${filename}' '${url}'
    """

    stub:
    """
    touch '${filename}'
    """
}

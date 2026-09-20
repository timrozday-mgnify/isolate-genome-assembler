// `plassembler download` fetches PLSDB and builds the mash sketches Plassembler needs.
// Part of `--prepare_databases`, so the pipeline proper only reads --plassembler_db.
process PLASSEMBLER_DOWNLOAD {
    tag 'plassembler_db'
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/plassembler:1.8.5--pyhdfd78af_0'
        : 'quay.io/biocontainers/plassembler:1.8.5--pyhdfd78af_0'}"

    output:
    tuple val('plassembler_db'), path('plassembler_db'), emit: database
    tuple val("${task.process}"), val('plassembler'), eval("plassembler --version | sed 's/^.*version //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plassembler download -d plassembler_db
    """

    stub:
    """
    mkdir -p plassembler_db
    touch plassembler_db/plsdb.msh
    """
}

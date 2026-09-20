// Build the human sylph database from CHM13 v2.0 and GRCh38. Human contamination is
// usually low-abundance, so it is screened with `sylph query` against this small database
// rather than with `sylph profile`.
process SYLPH_SKETCH_GENOMES {
    tag "${name}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/sylph:0.9.0--ha6fb395_0'
        : 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'}"

    input:
    tuple val(name), path(genomes)

    output:
    tuple val(name), path("${name}.syldb"), emit: database
    tuple val("${task.process}"), val('sylph'), eval('sylph -V | sed "s/sylph //g"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sylph sketch -g ${genomes} -t ${task.cpus} ${args} -o ${name}
    """

    stub:
    """
    touch ${name}.syldb
    """
}

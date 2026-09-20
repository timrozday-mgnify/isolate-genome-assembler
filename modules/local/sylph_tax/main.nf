// Attach GTDB taxonomy to a sylph profile. sylph reports genome accessions only, so the
// species names the gates and the report use come from here.
process SYLPH_TAX {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/sylph-tax:1.9.1--pyhdfd78af_0'
        : 'quay.io/biocontainers/sylph-tax:1.9.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(profile)
    path taxonomy

    output:
    tuple val(meta), path("${meta.id}.sylph_tax.tsv"), emit: taxonomy
    tuple val("${task.process}"), val('sylph-tax'), eval("sylph-tax --version | sed 's/^.* //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sylph-tax taxprof ${profile} \\
        -t ${taxonomy} \\
        ${args} \\
        -o ./
    mv *.sylphmpa ${meta.id}.sylph_tax.tsv
    """

    stub:
    """
    printf 'clade_name\\trelative_abundance\\tsequence_abundance\\nd__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia coli\\t100.0\\t100.0\\n' > ${meta.id}.sylph_tax.tsv
    """
}

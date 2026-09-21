// Myloasm's polymorphic-k-mer string graph: fast, and a different algorithm from the
// default set. Off by default; flags follow Autocycler v0.7.0's helper.rs.
process MYLOASM {
    tag "${meta.id}_${subset}"
    label 'process_medium'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/myloasm:0.7.0--hdcadc20_0'
        : 'quay.io/biocontainers/myloasm:0.7.0--hdcadc20_0'}"

    input:
    tuple val(meta), val(subset), path(reads)

    output:
    tuple val(meta), val('myloasm'), val(subset), path('out/*'), emit: assembly
    tuple val("${task.process}"), val('myloasm'), eval("myloasm --version | sed 's/^myloasm //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    myloasm --output-dir myloasm_out ${reads} --threads ${task.cpus} --hifi ${args}
    mkdir -p out
    cp myloasm_out/assembly_primary.fa out/
    """

    stub:
    """
    mkdir -p out
    printf '>u1ctg_len-4_circular-yes_depth-30-30-30_duplicated-no mult=1.00\\nACGT\\n' > out/assembly_primary.fa
    """
}

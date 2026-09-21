// hifiasm's string graph, built for HiFi reads. Flags follow Autocycler v0.7.0's helper.rs.
process HIFIASM {
    tag "${meta.id}_${subset}"
    label 'process_high'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/hifiasm:0.25.0--h5ca1c30_0'
        : 'quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0'}"

    input:
    tuple val(meta), val(subset), path(reads)

    output:
    tuple val(meta), val('hifiasm'), val(subset), path('out/*'), emit: assembly
    tuple val("${task.process}"), val('hifiasm'), eval("hifiasm --version"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p out
    hifiasm -t ${task.cpus} -o out/hifiasm -l 0 -f 0 ${args} ${reads}
    # hifiasm writes many graph files; the primary contig graph is the assembly.
    find out -mindepth 1 -maxdepth 1 ! -name 'hifiasm.bp.p_ctg.gfa' -exec rm -rf {} +
    """

    stub:
    """
    mkdir -p out
    printf 'S\\tptg000001c\\tACGT\\tdp:f:30\\n' > out/hifiasm.bp.p_ctg.gfa
    """
}

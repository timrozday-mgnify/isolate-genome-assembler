// Flye's repeat graph, one of the seven input assemblies. Commands follow Autocycler
// v0.7.0's src/helper.rs; the native output goes to out/ and NORMALISE_HEADERS converts it.
process FLYE {
    tag "${meta.id}_${subset}"
    label 'process_high'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/flye:2.9.6--py312h734f728_1'
        : 'quay.io/biocontainers/flye:2.9.6--py312h734f728_1'}"

    input:
    tuple val(meta), val(subset), path(reads)

    output:
    tuple val(meta), val('flye'), val(subset), path('out/*'), emit: assembly
    tuple val("${task.process}"), val('flye'), eval("flye --version"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // ext.allow_no_assembly: reads too short or too few to overlap give an empty assembly_info.txt, not a failure.
    def on_no_assembly = task.ext.allow_no_assembly
        ? "|| { status=\$?; grep -qE 'No disjointigs were assembled|No reads above minimum length threshold' out/flye.log || exit \$status; touch out/assembly_info.txt; }"
        : ''
    """
    flye --pacbio-hifi ${reads} --threads ${task.cpus} --out-dir out ${args} ${on_no_assembly}
    # Only the files NORMALISE_HEADERS and the report need are kept.
    find out -mindepth 1 -maxdepth 1 ! -name assembly.fasta ! -name assembly_info.txt \\
        ! -name assembly_graph.gfa ! -name flye.log -exec rm -rf {} +
    """

    stub:
    """
    mkdir -p out
    printf '>contig_1\\nACGT\\n' > out/assembly.fasta
    printf '#seq_name\\tlength\\tcov.\\tcirc.\\n contig_1\\t4\\t30\\tY\\n' > out/assembly_info.txt
    """
}

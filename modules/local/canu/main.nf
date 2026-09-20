// Canu's overlap-consensus assembly. It is the slowest of the seven, hence process_long;
// -fast and useGrid=false follow Autocycler v0.7.0's helper.rs.
process CANU {
    tag "${meta.id}_${subset}"
    label 'process_high'
    label 'process_long'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/canu:2.3--h636b4d1_3'
        : 'quay.io/biocontainers/canu:2.3--h636b4d1_3'}"

    input:
    tuple val(meta), val(subset), path(reads), val(genome_size)

    output:
    tuple val(meta), val('canu'), val(subset), path('out/*'), emit: assembly
    tuple val("${task.process}"), val('canu'), eval("canu -version | sed 's/^Canu //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    canu -p canu -d out \\
        -fast \\
        genomeSize=${genome_size} \\
        useGrid=false \\
        maxThreads=${task.cpus} \\
        ${args} \\
        -pacbio-hifi ${reads}

    # The tigInfo file carries the depths NORMALISE_HEADERS puts in the headers.
    find out -mindepth 1 -maxdepth 1 ! -name 'canu.contigs.fasta' \\
        ! -name 'canu.contigs.layout.tigInfo' ! -name 'canu.report' -exec rm -rf {} +
    """

    stub:
    """
    mkdir -p out
    printf '>tig00000001 len=4 suggestCircular=yes trim=0-4\\nACGT\\n' > out/canu.contigs.fasta
    printf '#tigID\\tlength\\tcoverage\\n1\\t4\\t30\\n' > out/canu.contigs.layout.tigInfo
    """
}

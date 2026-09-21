// Plassembler as an Autocycler input: it finds small plasmids the whole-genome assemblers
// miss, and NORMALISE_HEADERS gives its circular contigs double cluster weight so a
// plasmid only it saw still passes cluster QC.
process PLASSEMBLER {
    tag "${meta.id}_${subset}"
    label 'process_medium'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/plassembler:1.8.5--pyhdfd78af_0'
        : 'quay.io/biocontainers/plassembler:1.8.5--pyhdfd78af_0'}"

    input:
    tuple val(meta), val(subset), path(reads)
    path plassembler_db

    output:
    tuple val(meta), val('plassembler'), val(subset), path('out/*'), emit: assembly
    tuple val(meta), path('out/plassembler_plasmids.fasta'), emit: plasmids
    tuple val(meta), path('out/plassembler_summary.tsv'), emit: summary, optional: true
    tuple val("${task.process}"), val('plassembler'), eval("plassembler --version | sed 's/^.*version //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Plassembler exits if `canu --version` fails, even though canu is never run here. Under
    // Singularity on some hosts the image's relocatable Perl resolves @INC to ../lib/perl5
    // relative to the working directory and cannot load even strict.pm, so the paths are
    // given explicitly.
    """
    export PERL5LIB=/usr/local/lib/perl5/5.32/site_perl:/usr/local/lib/perl5/site_perl:/usr/local/lib/perl5/5.32/vendor_perl:/usr/local/lib/perl5/vendor_perl:/usr/local/lib/perl5/5.32/core_perl:/usr/local/lib/perl5/core_perl

    status=0
    plassembler long \\
        -d ${plassembler_db} \\
        -l ${reads} \\
        -o out \\
        -t ${task.cpus} \\
        --force \\
        --skip_qc \\
        --pacbio_model pacbio-hifi \\
        ${args} 2> plassembler.log || status=\$?
    cat plassembler.log >&2

    # Plassembler exits 1 when its Flye run yields no chromosome-length contig, which shallow
    # reads make routine. It cannot separate plasmids from chromosome then, so this is
    # reported as no plasmids rather than a failure.
    if [ "\$status" -ne 0 ]; then
        grep -q 'No chromosome was identified' plassembler.log || exit "\$status"
        rm -rf out
        mkdir out
    fi

    find out -mindepth 1 -maxdepth 1 ! -name 'plassembler_plasmids.fasta' \\
        ! -name 'plassembler_plasmids.gfa' ! -name 'plassembler_summary.tsv' -exec rm -rf {} +
    # An isolate with no plasmids is a normal result, not a failure.
    touch out/plassembler_plasmids.fasta
    """

    stub:
    """
    mkdir -p out
    printf '>1 len=4 circular=true\\nACGT\\n' > out/plassembler_plasmids.fasta
    printf 'contig\\tlength\\tcopy_number_long\\tPLSDB_hit\\n1\\t4\\t1.0\\tNone\\n' > out/plassembler_summary.tsv
    """
}

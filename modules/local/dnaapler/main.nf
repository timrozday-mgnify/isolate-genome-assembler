// Rotate each contig to a consistent start: the chromosome at dnaA, plasmids at repA,
// phages at terL. `--autocomplete mystery` gives a contig with no gene hit a reproducible
// arbitrary start rather than leaving it unrotated, so assemblies stay comparable.
process DNAAPLER {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/dnaapler:1.4.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/dnaapler:1.4.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}_reoriented.fasta"), emit: rotated
    tuple val(meta), path("${meta.id}_failed_to_reorient.fasta"), emit: unrotated, optional: true
    tuple val("${task.process}"), val('dnaapler'), eval("dnaapler --version | sed 's/^.*version //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--autocomplete mystery'
    """
    # Every --autocomplete mode exits dnaapler outright on a no-hit contig with fewer than
    # 4 CDS (mystery/largest; nearest needs 2), so list those small contigs for --ignore.
    # They are kept in the output as-is. Counted with the same pyrodigal call dnaapler uses.
    # ponytail: a small contig that does have a gene hit is ignored too; rotation of a <4 CDS
    # contig is cosmetic.
    python3 - <<'EOF' > ignore.txt
    import pyrodigal
    from Bio import SeqIO
    finder = pyrodigal.GeneFinder(meta=True)
    for record in SeqIO.parse("${assembly}", "fasta"):
        if len(finder.find_genes(str(record.seq))) < 4:
            print(record.id)
    EOF

    dnaapler all \\
        --input ${assembly} \\
        --output dnaapler \\
        --prefix ${meta.id} \\
        --threads ${task.cpus} \\
        --ignore ignore.txt \\
        --force \\
        ${args}

    # dnaapler writes nothing when every contig was already oriented; the classifier still
    # needs a file to read.
    touch dnaapler/${meta.id}_reoriented.fasta
    mv dnaapler/${meta.id}_reoriented.fasta ${meta.id}_reoriented.fasta
    if [ -s dnaapler/${meta.id}_failed_to_reorient.fasta ]; then
        mv dnaapler/${meta.id}_failed_to_reorient.fasta ${meta.id}_failed_to_reorient.fasta
    fi
    """

    stub:
    """
    printf '>contig_1 length=4 depth=30 circular=true\\nACGT\\n' > ${meta.id}_reoriented.fasta
    """
}

// An assembly the samplesheet supplied, put into the shape the finishing stage expects:
// decompressed, with blank lines and Windows line endings gone.
//
// Nothing is added. In particular no `circular=` tag is invented: a submitted FASTA does
// not say whether a contig is closed, and classify_replicons.py reports that as unknown
// rather than letting it read as linear.
process PREPARE_ASSEMBLY {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), path(assembly, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}.input_assembly.fasta"), emit: assembly

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gzip -cdf ${assembly} | tr -d '\\r' | grep -v '^\$' > ${meta.id}.input_assembly.fasta

    if ! grep -q '^>' ${meta.id}.input_assembly.fasta; then
        echo "ERROR: ${assembly} holds no FASTA records" >&2
        exit 1
    fi
    """

    stub:
    """
    printf '>contig_1 length=4\\nACGT\\n' > ${meta.id}.input_assembly.fasta
    """
}

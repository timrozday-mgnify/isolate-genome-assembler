// One assembler's native output becomes an Autocycler input assembly: contigs renamed
// <assembler>_<subset>_<n> so a consensus contig can be traced back, with the circularity
// and depth tags each tool reports. The per-assembler rules live in bin/normalise_headers.py.
process NORMALISE_HEADERS {
    tag "${meta.id}_${assembler}_${subset}"
    label 'process_single'
    label 'assembler'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), val(assembler), val(subset), path(assembly, stageAs: 'native/*')

    output:
    tuple val(meta), path("${assembler}_${subset}.fasta"), emit: assembly
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    normalise_headers.py \\
        --assembler ${assembler} \\
        --subset ${subset} \\
        --input-dir native \\
        --output ${assembler}_${subset}.fasta
    """

    stub:
    def fail = "${assembler}_${subset}" in params.stub_fail_assemblies.tokenize(',')
    """
    ${fail ? 'exit 1' : ''}
    printf '>${assembler}_${subset}_1 length=4 circular=true\\nACGT\\n' > ${assembler}_${subset}.fasta
    """
}

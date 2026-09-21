// Allele depths at every site where reads disagree with the assembly (check D). There is
// deliberately no `bcftools call`: a haploid caller reports only the majority allele and
// would hide the minority alleles that point to a mixed strain or a collapsed repeat.
process VARIANT_PILEUP {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.23.1--hb2cee57_0'
        : 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'}"

    input:
    tuple val(meta), path(bam), path(bai), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.pileup.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools mpileup \\
        --fasta-ref ${assembly} \\
        --annotate FORMAT/AD,FORMAT/DP \\
        --threads ${task.cpus} \\
        -Ou \\
        ${args} \\
        ${bam} \\
        | bcftools view -i 'FMT/AD[0:1] >= 2' -Oz -o ${meta.id}.pileup.vcf.gz
    """

    stub:
    """
    echo '##fileformat=VCFv4.2' | bgzip > ${meta.id}.pileup.vcf.gz
    """
}

// Which of the sample x subset x assembler jobs produced an input assembly. Assembler
// processes are allowed to fail (one crashing on one subset must not sink the sample), so
// the missing rows here are the record of that, alongside the Nextflow trace.
process ASSEMBLY_ATTEMPTS {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(meta), val(expected), path(assemblies)

    output:
    tuple val(meta), path("${meta.id}.assembly_attempts.tsv"), emit: attempts

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    printf 'sample\\tassembler\\tsubset\\tstatus\\tcontigs\\n' > ${meta.id}.assembly_attempts.tsv
    for attempt in ${expected.join(' ')}; do
        assembler=\${attempt%_*}
        subset=\${attempt##*_}
        if [ -s "\${attempt}.fasta" ]; then
            contigs=\$(grep -c '^>' "\${attempt}.fasta" || true)
            status=ok
        else
            contigs=0
            status=failed
        fi
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' \\
            '${meta.id}' "\$assembler" "\$subset" "\$status" "\$contigs" \\
            >> ${meta.id}.assembly_attempts.tsv
    done
    """

    stub:
    """
    printf 'sample\\tassembler\\tsubset\\tstatus\\tcontigs\\n' > ${meta.id}.assembly_attempts.tsv
    for attempt in ${expected.join(' ')}; do
        assembler=\${attempt%_*}
        subset=\${attempt##*_}
        if [ -s "\${attempt}.fasta" ]; then
            contigs=\$(grep -c '^>' "\${attempt}.fasta" || true)
            status=ok
        else
            contigs=0
            status=failed
        fi
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' \\
            '${meta.id}' "\$assembler" "\$subset" "\$status" "\$contigs" \\
            >> ${meta.id}.assembly_attempts.tsv
    done
    """
}

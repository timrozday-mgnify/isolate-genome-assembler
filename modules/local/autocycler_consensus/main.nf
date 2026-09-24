// Autocycler steps 3-7 in one process: compress, cluster, then trim/resolve per QC-pass
// cluster and combine. They all read and write one shared autocycler_out/ directory and
// finish in minutes, so splitting them would only mean staging that directory around.
//
// The process is deliberately failure-tolerant. `cluster` refuses to run past
// --max_contigs and `resolve` can give up on a tangled cluster; either way the sample
// still gets an assembly, because SELECT_ASSEMBLY falls back to the full-read Flye
// assembly whenever combine did not report a fully resolved consensus. Autocycler's own
// output and log are published either way.
process AUTOCYCLER_CONSENSUS {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/autocycler:0.7.0--h79ce301_0'
        : 'quay.io/biocontainers/autocycler:0.7.0--h79ce301_0'}"

    input:
    tuple val(meta), path(assemblies, stageAs: 'assemblies/*'), path(reads), path(autocycler_in, stageAs: 'autocycler_in')

    output:
    tuple val(meta), path('autocycler_out'), emit: autocycler_dir
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: consensus
    tuple val(meta), path("${meta.id}.consensus.yaml"), emit: metrics
    tuple val(meta), path("${meta.id}.autocycler_table.tsv"), emit: table
    tuple val(meta), path('dotplots/*.png'), emit: dotplots, optional: true
    tuple val(meta), path('autocycler_out/consensus_assembly.gfa'), emit: gfa, optional: true
    tuple val("${task.process}"), val('autocycler'), eval("autocycler --version | sed 's/^.*r //; s/^autocycler //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cluster_args = meta.autocycler_cluster_args ?: ''
    """
    set +e
    mkdir -p autocycler_out
    {
        if [ -d autocycler_in ]; then
            # Re-entry: a previous run's directory is reclustered, so the expensive
            # compress step is not repeated.
            cp -rL autocycler_in autocycler_out
            rm -rf autocycler_out/clustering
        else
            autocycler compress -i assemblies -a autocycler_out -t ${task.cpus} ${args}
        fi

        autocycler cluster -a autocycler_out ${cluster_args}

        for cluster in autocycler_out/clustering/qc_pass/cluster_*; do
            autocycler trim -c "\$cluster" -t ${task.cpus}
            autocycler resolve -c "\$cluster"
        done

        autocycler combine \\
            -a autocycler_out \\
            -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa \\
            --reads ${reads} \\
            -t ${task.cpus}
    } 2>&1 | tee ${meta.id}.autocycler.log
    set -e

    mkdir -p dotplots
    for cluster in autocycler_out/clustering/qc_pass/cluster_*; do
        [ -f "\$cluster/5_final.gfa" ] || continue
        autocycler dotplot -i "\$cluster/5_final.gfa" \\
            -o "dotplots/${meta.id}_\$(basename "\$cluster").png" \\
            || echo "dotplot failed for \$cluster" >&2
    done

    if [ -f autocycler_out/consensus_assembly.yaml ]; then
        cp autocycler_out/consensus_assembly.yaml ${meta.id}.consensus.yaml
        cp autocycler_out/consensus_assembly.fasta ${meta.id}.consensus.fasta
        # `table -a` writes the row only, so the header comes from a bare `table`.
        # Without it collect_metrics.py reads the row as the header and the sample
        # vanishes from the report.
        autocycler table > ${meta.id}.autocycler_table.tsv
        autocycler table -a autocycler_out -n ${meta.id} >> ${meta.id}.autocycler_table.tsv
    else
        # Autocycler stopped before combine. Say so in its own vocabulary, so
        # select_assembly.py needs no special case for a crash.
        echo 'consensus_assembly_fully_resolved: false' > ${meta.id}.consensus.yaml
        touch ${meta.id}.consensus.fasta
        autocycler table > ${meta.id}.autocycler_table.tsv
    fi
    mv ${meta.id}.autocycler.log autocycler_out/
    """

    stub:
    """
    mkdir -p autocycler_out/clustering/qc_pass/cluster_001
    printf '>cluster_1 length=4 circular=true\\nACGT\\n' > ${meta.id}.consensus.fasta
    cp ${meta.id}.consensus.fasta autocycler_out/consensus_assembly.fasta
    printf 'H\\tVN:Z:1.0\\nS\\t1\\tACGT\\n' > autocycler_out/consensus_assembly.gfa
    echo 'consensus_assembly_fully_resolved: ${params.stub_fully_resolved}' > ${meta.id}.consensus.yaml
    cp ${meta.id}.consensus.yaml autocycler_out/consensus_assembly.yaml
    printf 'name\\tconsensus_assembly_fully_resolved\\n${meta.id}\\t${params.stub_fully_resolved}\\n' > ${meta.id}.autocycler_table.tsv
    mkdir -p dotplots
    printf '\\211PNG\\r\\n\\032\\n' > dotplots/${meta.id}_cluster_001.png
    """
}

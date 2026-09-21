// Map the reads back to the finished assembly. Everything in check B reads this BAM, and
// the reads that do not map go on to their own quick assembly: a circular contig among
// them is a replicon the assembly may be missing.
//
// Circular contigs are mapped as they are. A read spanning the origin is split into a
// primary and a supplementary alignment there, which mosdepth counts in full and the
// clipping check ignores, so no rotated or extended copy of the reference is needed.
process MAP_READS {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data'
        : 'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce'}"

    input:
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.mapping.tsv"), emit: stats
    tuple val(meta), path("${meta.id}.unmapped.fastq.gz"), emit: unmapped
    tuple val("${task.process}"), val('minimap2'), eval('minimap2 --version'), topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d; s/.* //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 -t ${task.cpus} -ax map-hifi --secondary=no ${args} ${assembly} ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.bam
    samtools index ${meta.id}.bam
    samtools fastq -f 4 -0 ${meta.id}.unmapped.fastq.gz ${meta.id}.bam

    # One primary record per read (-F 0x900), so these are read and base fractions.
    samtools view -F 0x900 ${meta.id}.bam | awk -v sample=${meta.id} '
        BEGIN { OFS = "\\t" }
        {
            reads++; bases += length(\$10)
            if (int(\$2 / 4) % 2) { unmapped++; unmapped_bases += length(\$10) }
        }
        END {
            print "sample", "reads", "unmapped_reads", "unmapped_read_fraction", "bases", "unmapped_bases", "unmapped_base_fraction"
            print sample, reads + 0, unmapped + 0, (reads ? unmapped / reads : 0), bases + 0, unmapped_bases + 0, (bases ? unmapped_bases / bases : 0)
        }' > ${meta.id}.mapping.tsv
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai
    printf '@read_001\\nACGT\\n+\\n####\\n' | bgzip > ${meta.id}.unmapped.fastq.gz
    printf 'sample\\treads\\tunmapped_reads\\tunmapped_read_fraction\\tbases\\tunmapped_bases\\tunmapped_base_fraction\\n${meta.id}\\t1\\t0\\t0\\t4\\t0\\t0\\n' > ${meta.id}.mapping.tsv
    """
}

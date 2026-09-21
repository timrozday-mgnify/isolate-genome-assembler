include { MAP_READS                   } from '../../modules/local/map_reads'
include { MOSDEPTH                    } from '../../modules/nf-core/mosdepth'
include { COVERAGE_REGIONS            } from '../../modules/local/coverage_regions'
include { CLIPPING_PILEUPS            } from '../../modules/local/clipping_pileups'
include { INSPECTOR                   } from '../../modules/local/inspector'
include { FLYE as FLYE_UNMAPPED       } from '../../modules/local/flye'
include { VARIANT_PILEUP              } from '../../modules/local/variant_pileup'
include { VARIANT_SCAN                } from '../../modules/local/variant_scan'
include { MERYL_COUNT                 } from '../../modules/nf-core/meryl/count'
include { MERQURY_MERQURY             } from '../../modules/nf-core/merqury/merqury'
include { BAKTA_BAKTA                 } from '../../modules/nf-core/bakta/bakta'
include { DIAMOND_BLASTP              } from '../../modules/nf-core/diamond/blastp'
include { GENE_CHECKS                 } from '../../modules/local/gene_checks'
include { CHECKM2_PREDICT             } from '../../modules/nf-core/checkm2/predict'
include { BUSCO_BUSCO                 } from '../../modules/nf-core/busco/busco'
include { GTDBTK_CLASSIFYWF           } from '../../modules/nf-core/gtdbtk/classifywf'
include { SKANI_DIST as SKANI_REFERENCE } from '../../modules/nf-core/skani/dist'
include { MINIMAP2_ALIGN as MINIMAP2_REFERENCE } from '../../modules/nf-core/minimap2/align'
include { BANDAGE_IMAGE               } from '../../modules/nf-core/bandage/image'
include { ASSEMBLER_CONTRIBUTION      } from '../../modules/local/assembler_contribution'

// `main.nf` requires these params whenever there are samples; returning an empty path
// keeps a sample-free `-preview` from tripping over a null.
def database(param) {
    param ? file(param, checkIfExists: true) : []
}

// Merqury's recommended k for a genome of this size (its best_k.sh): the shortest k at
// which a random k-mer is unlikely to recur, with a 0.1% collision rate.
def merylK(genomeSize) {
    def size = genomeSize?.toString()?.isDouble() ? genomeSize as double : 5_000_000
    Math.max(15, Math.ceil(Math.log(size * (1 - 0.001) / 0.001) / Math.log(4)) as int)
}

// Stage 6. Every check writes a small table that QC_GATES and the report read; none of
// them decides pass or fail itself.
workflow CHECKS {
    take:
    ch_assembly // tuple: [meta, final assembly FASTA]
    ch_contigs // tuple: [meta, <id>.contigs.tsv]
    ch_reads // tuple: [meta, screened HiFi FASTQ]
    ch_read_qc // tuple: [meta, read_qc.tsv], for the genome size
    ch_autocycler_dir // tuple: [meta, autocycler_out/]
    ch_consensus_gfa // tuple: [meta, consensus_assembly.gfa], when combine ran

    main:
    // --- B. Read support ---
    MAP_READS(ch_assembly.join(ch_reads))
    MOSDEPTH(MAP_READS.out.bam.map { meta, bam, bai -> [meta, bam, bai, []] }, [[:], []], false)
    COVERAGE_REGIONS(MOSDEPTH.out.regions_bed.join(ch_contigs))
    CLIPPING_PILEUPS(MAP_READS.out.bam)
    INSPECTOR(ch_assembly.join(ch_reads))

    // A circular contig among the reads that did not map is a replicon the assembly may
    // be missing. Flye cannot assemble nothing, so a sample whose reads all mapped skips it.
    FLYE_UNMAPPED(
        MAP_READS.out.unmapped
            .filter { _meta, fastq -> fastq.countFastq() > 0 }
            .map { meta, fastq -> [meta, 'unmapped', fastq] }
    )
    ch_unmapped_info = FLYE_UNMAPPED.out.assembly.map { meta, _assembler, _subset, files ->
        [meta, files.find { it.name == 'assembly_info.txt' }]
    }

    // --- D. Per-base accuracy ---
    VARIANT_PILEUP(MAP_READS.out.bam.join(ch_assembly))
    VARIANT_SCAN(VARIANT_PILEUP.out.vcf.join(ch_assembly))

    ch_meryl_in = ch_reads
        .join(ch_read_qc.splitCsv(header: true, sep: '\t').map { meta, row -> [meta, row.genome_size_used] })
        .multiMap { meta, reads, genome_size ->
            reads: [meta, reads]
            k: merylK(genome_size)
        }
    MERYL_COUNT(ch_meryl_in.reads, ch_meryl_in.k)
    MERQURY_MERQURY(MERYL_COUNT.out.meryl_db.join(ch_assembly))

    // --- E. Gene-level integrity ---
    BAKTA_BAKTA(ch_assembly, database(params.bakta_db), [], [], [], [])
    DIAMOND_BLASTP(
        BAKTA_BAKTA.out.faa,
        [[id: 'ideel_db'], database(params.ideel_db)],
        6,
        'qseqid sseqid qlen slen pident bitscore',
    )
    // DIAMOND_BLASTP adds the database name to meta, which would stop every later join.
    ch_ideel_hits = DIAMOND_BLASTP.out.txt.map { meta, hits -> [meta - [db: 'ideel_db'], hits] }
    GENE_CHECKS(
        BAKTA_BAKTA.out.faa
            .join(ch_ideel_hits)
            .join(BAKTA_BAKTA.out.gff)
            .join(MOSDEPTH.out.regions_bed)
    )

    // --- F. Completeness, contamination and identity ---
    CHECKM2_PREDICT(ch_assembly, [[id: 'checkm2_db'], database(params.checkm2_db)])
    BUSCO_BUSCO(ch_assembly, 'genome', params.busco_lineage, database(params.busco_db), [], true)
    GTDBTK_CLASSIFYWF(ch_assembly, ['gtdbtk_db', database(params.gtdbtk_db)], false)

    // skani against the samplesheet reference, for the samples that gave one. multiMap
    // keeps query and reference in lockstep across samples.
    ch_reference_in = ch_assembly
        .filter { meta, _assembly -> meta.reference }
        .multiMap { meta, assembly ->
            query: [meta, assembly]
            reference: [meta, meta.reference]
        }
    SKANI_REFERENCE(ch_reference_in.query, ch_reference_in.reference)
    // The alignment blocks the report draws as a dotplot against the reference.
    MINIMAP2_REFERENCE(ch_reference_in.query, ch_reference_in.reference, false, '', false, false)

    // --- A and C. Graph structure and consensus consistency ---
    BANDAGE_IMAGE(ch_consensus_gfa)
    ASSEMBLER_CONTRIBUTION(ch_autocycler_dir)

    // Every measurement QC_GATES reads, as [meta, qc_gates.py option, file].
    metrics = channel.empty().mix(
        MAP_READS.out.stats.map { meta, f -> [meta, 'mapping', f] },
        COVERAGE_REGIONS.out.regions.map { meta, f -> [meta, 'coverage-regions', f] },
        COVERAGE_REGIONS.out.depth.map { meta, f -> [meta, 'contig-depth', f] },
        CLIPPING_PILEUPS.out.pileups.map { meta, f -> [meta, 'clipping', f] },
        INSPECTOR.out.summary.map { meta, f -> [meta, 'inspector', f] },
        ch_unmapped_info.map { meta, f -> [meta, 'unmapped-assembly', f] },
        VARIANT_SCAN.out.summary.map { meta, f -> [meta, 'variants', f] },
        MERQURY_MERQURY.out.assembly_qv.map { meta, f -> [meta, 'merqury', f] },
        GENE_CHECKS.out.summary.map { meta, f -> [meta, 'ideel', f] },
        CHECKM2_PREDICT.out.checkm2_tsv.map { meta, f -> [meta, 'checkm2', f] },
        BUSCO_BUSCO.out.batch_summary.map { meta, f -> [meta, 'busco', f] },
        // One summary per domain GTDB-Tk found markers for: a single path or a list.
        GTDBTK_CLASSIFYWF.out.summary.flatMap { meta, fs -> (fs instanceof List ? fs : [fs]).collect { f -> [meta, 'gtdbtk', f] } },
    )

    // What only the report reads, in the same [meta, kind, file] shape; see
    // bin/collect_metrics.py for how each kind is read.
    report = channel.empty().mix(
        MOSDEPTH.out.regions_bed.map { meta, f -> [meta, 'coverage', f] },
        VARIANT_SCAN.out.variants.map { meta, f -> [meta, 'variant-sites', f] },
        MERQURY_MERQURY.out.stats.map { meta, f -> [meta, 'merqury-completeness', f] },
        MERQURY_MERQURY.out.spectra_cn_fl_png.map { meta, f -> [meta, 'image-spectra-cn', f] },
        BAKTA_BAKTA.out.txt.map { meta, f -> [meta, 'bakta', f] },
        GENE_CHECKS.out.ratios.map { meta, f -> [meta, 'ideel-ratios', f] },
        GENE_CHECKS.out.rrna.map { meta, f -> [meta, 'rrna-depth', f] },
        SKANI_REFERENCE.out.dist.map { meta, f -> [meta, 'reference-skani', f] },
        MINIMAP2_REFERENCE.out.paf.map { meta, f -> [meta, 'dotplot', f] },
        BANDAGE_IMAGE.out.png.map { meta, f -> [meta, 'image-bandage', f] },
        ASSEMBLER_CONTRIBUTION.out.contribution.map { meta, f -> [meta, 'assembler-contribution', f] },
    )

    emit:
    metrics
    report
    bam = MAP_READS.out.bam
}

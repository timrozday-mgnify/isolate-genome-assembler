include { PREPARE_READS          } from '../../modules/local/prepare_reads'
include { HIFIADAPTERFILT        } from '../../modules/local/hifiadapterfilt'
include { READ_PROFILE           } from '../../modules/local/read_profile'
include { AUTOCYCLER_GENOME_SIZE } from '../../modules/local/autocycler_genome_size'
include { KMC_HISTOGRAM          } from '../../modules/local/kmc_histogram'
include { GENOMESCOPE2           } from '../../modules/local/genomescope2'
include { READ_QC_SUMMARY        } from '../../modules/local/read_qc_summary'
include { SEQKIT_STATS           } from '../../modules/nf-core/seqkit/stats'
include { NANOPLOT               } from '../../modules/nf-core/nanoplot'

// Stage 1. Normalise the inputs, measure them, and hand the screened read set on. Every
// measurement ends up as a column of read_qc.tsv; none of them is turned into a status here.
workflow READ_QC {
    take:
    ch_samples // tuple: [meta, one or more HiFi read files]

    main:
    PREPARE_READS(ch_samples)
    HIFIADAPTERFILT(PREPARE_READS.out.reads)

    // The adapter-filtered set is what the rest of the pipeline sees, unless removal is off.
    ch_reads = params.remove_adapter_reads ? HIFIADAPTERFILT.out.reads : PREPARE_READS.out.reads

    SEQKIT_STATS(ch_reads)
    NANOPLOT(ch_reads)
    READ_PROFILE(ch_reads)
    AUTOCYCLER_GENOME_SIZE(ch_reads)
    KMC_HISTOGRAM(ch_reads)
    GENOMESCOPE2(KMC_HISTOGRAM.out.histogram)

    READ_QC_SUMMARY(
        PREPARE_READS.out.stats
            .join(SEQKIT_STATS.out.stats)
            .join(READ_PROFILE.out.gc_hist)
            .join(READ_PROFILE.out.duplicates)
            .join(HIFIADAPTERFILT.out.stats)
            .join(AUTOCYCLER_GENOME_SIZE.out.genome_size)
            .join(GENOMESCOPE2.out.summary)
    )

    emit:
    reads = ch_reads
    summary = READ_QC_SUMMARY.out.summary
    nanoplot = NANOPLOT.out.txt
    genomescope = GENOMESCOPE2.out.results
}

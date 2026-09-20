include { SYLPH_PROFILE         } from '../../modules/nf-core/sylph/profile'
include { SYLPH_QUERY           } from '../../modules/nf-core/sylph/query'
include { MINIMAP2_ALIGN        } from '../../modules/nf-core/minimap2/align'
include { SYLPH_TAX             } from '../../modules/local/sylph_tax'
include { HUMAN_FRACTION        } from '../../modules/local/human_fraction'
include { CONTAMINATION_SUMMARY } from '../../modules/local/contamination_summary'

// `main.nf` requires these params whenever there are samples; returning an empty path
// keeps a sample-free `-preview` from tripping over a null.
def database(param) {
    param ? file(param, checkIfExists: true) : []
}

// Stage 2. Species profile against GTDB, plus two independent human measures: sylph's
// containment query (sensitive at low abundance) and the mapped read fraction.
workflow CONTAMINATION {
    take:
    ch_reads // tuple: [meta, normalised HiFi FASTQ]

    main:
    // sylph's read flags are chosen from meta.single_end; HiFi is always single-end.
    ch_single = ch_reads.map { meta, reads -> [meta + [single_end: true], reads] }

    SYLPH_PROFILE(ch_single, database(params.sylph_gtdb_db))
    SYLPH_TAX(SYLPH_PROFILE.out.profile_out, database(params.sylph_gtdb_taxonomy))
    SYLPH_QUERY(ch_single, database(params.sylph_human_db))

    MINIMAP2_ALIGN(
        ch_reads,
        [[id: 'human_reference'], database(params.human_reference)],
        false, // bam_format
        '', // bam_index_extension
        false, // cigar_paf_format
        false, // cigar_bam
    )
    HUMAN_FRACTION(MINIMAP2_ALIGN.out.paf.join(ch_reads))

    ch_sylph = SYLPH_PROFILE.out.profile_out
        .join(SYLPH_TAX.out.taxonomy)
        .join(SYLPH_QUERY.out.query_out)
        .map { meta, profile, taxonomy, query -> [meta - [single_end: true], profile, taxonomy, query] }

    CONTAMINATION_SUMMARY(ch_sylph.join(HUMAN_FRACTION.out.fraction))

    emit:
    summary = CONTAMINATION_SUMMARY.out.summary
    remove_human = CONTAMINATION_SUMMARY.out.remove_human
    human_read_ids = HUMAN_FRACTION.out.read_ids
}

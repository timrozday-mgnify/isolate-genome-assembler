include { CONTIG_ENDS                    } from '../../modules/local/contig_ends'
include { DNAAPLER                       } from '../../modules/local/dnaapler'
include { CLASSIFY_REPLICONS             } from '../../modules/local/classify_replicons'
include { PLASSEMBLER as PLASSEMBLER_FULL } from '../../modules/local/plassembler'
include { PLASMID_AUDIT                  } from '../../modules/local/plasmid_audit'
include { SKANI_DIST                     } from '../../modules/nf-core/skani/dist'

// Stage 5. Turn the selected assembly into named, rotated replicons, and check it against
// an independent view of the sample's plasmids.
//
// The end-overlap check runs before rotation on purpose: rotating a contig moves a
// duplicated end into the middle of the sequence, where this check would no longer see it.
workflow FINISHING {
    take:
    ch_assembly // tuple: [meta, the assembly SELECT_ASSEMBLY chose]
    ch_reads // tuple: [meta, screened HiFi FASTQ], for the full-read plasmid audit

    main:
    CONTIG_ENDS(ch_assembly)
    DNAAPLER(ch_assembly)

    // dnaapler emits failed_to_reorient.fasta only when a contig had no gene hit at all,
    // which `--autocomplete mystery` makes rare.
    CLASSIFY_REPLICONS(
        DNAAPLER.out.rotated
            .join(DNAAPLER.out.unrotated, remainder: true)
            .map { meta, rotated, unrotated -> [meta, rotated, unrotated ?: []] }
            .join(CONTIG_ENDS.out.overlaps)
    )

    // Plassembler's second run is independent of the consensus: it sees the full read set,
    // so a plasmid it finds that the assembly lacks is real evidence of something missed.
    PLASSEMBLER_FULL(
        ch_reads.map { meta, reads -> [meta, 'full', reads] },
        params.plassembler_db ? file(params.plassembler_db, checkIfExists: true) : [],
    )

    // multiMap keeps the query and reference channels in lockstep, which a plain pair of
    // queue channels into SKANI_DIST would not guarantee once there is more than one sample.
    ch_skani_in = PLASSEMBLER_FULL.out.plasmids
        .join(CLASSIFY_REPLICONS.out.assembly)
        .multiMap { meta, plasmids, assembly ->
            query: [meta, plasmids]
            reference: [meta, assembly]
        }
    SKANI_DIST(ch_skani_in.query, ch_skani_in.reference)

    PLASMID_AUDIT(
        PLASSEMBLER_FULL.out.summary
            .join(SKANI_DIST.out.dist, remainder: true)
            .map { meta, summary, skani -> [meta, summary, skani ?: []] }
            .filter { _meta, _summary, skani -> skani }
    )

    emit:
    assembly = CLASSIFY_REPLICONS.out.assembly
    contigs = CLASSIFY_REPLICONS.out.table
    removed = CLASSIFY_REPLICONS.out.removed
    plasmid_audit = PLASMID_AUDIT.out.audit
}

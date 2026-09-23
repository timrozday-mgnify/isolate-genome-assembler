include { CONTIG_ENDS as CONTIG_ENDS_CANDIDATE } from '../../modules/local/contig_ends'
include { CIRCULARISE                          } from '../../modules/local/circularise'
include { SEQKIT_STATS as SEQKIT_STATS_CANDIDATE } from '../../modules/nf-core/seqkit/stats'
include { MERYL_COUNT                          } from '../../modules/nf-core/meryl/count'
include { MERQURY_MERQURY as MERQURY_CANDIDATE } from '../../modules/nf-core/merqury/merqury'
include { MAP_READS as MAP_READS_CANDIDATE     } from '../../modules/local/map_reads'
include { CLIPPING_PILEUPS as CLIPPING_PILEUPS_CANDIDATE } from '../../modules/local/clipping_pileups'
include { SCORE_ASSEMBLIES                     } from '../../modules/local/score_assemblies'
include { SELECT_ASSEMBLY                      } from '../../modules/local/select_assembly'

// Merqury's recommended k for a genome of this size (its best_k.sh): the shortest k at
// which a random k-mer is unlikely to recur, with a 0.1% collision rate.
def merylK(genomeSize) {
    def size = genomeSize?.toString()?.isDouble() ? genomeSize as double : 5_000_000
    Math.max(15, Math.ceil(Math.log(size * (1 - 0.001) / 0.001) / Math.log(4)) as int)
}

// Score every candidate assembly for a sample against the others, reference-free. The
// candidates are each assembler's full-read assembly plus the Autocycler consensus, and
// every signal comes from a check the pipeline already runs.
//
// Each candidate carries its own meta id, `<sample>.<assembler>`, because the checks name
// their outputs after it; `sample` keeps the original id for the regrouping at the end.
workflow SCORING {
    take:
    ch_candidates // tuple: [meta, assembler, assembly FASTA]
    ch_reads // tuple: [meta, screened HiFi FASTQ]
    ch_read_qc // tuple: [meta, read_qc.tsv], for the genome size
    ch_selection_in // tuple: [meta, consensus, consensus metrics, Flye fallback]

    main:
    ch_genome_size = ch_read_qc
        .splitCsv(header: true, sep: '\t')
        .map { meta, row -> [meta, row.genome_size_used] }

    // One meryl database per sample, built here and passed to CHECKS, which used to build
    // its own. Counting the reads once is the whole reason this moved.
    ch_meryl_in = ch_reads
        .join(ch_genome_size)
        .multiMap { meta, reads, genome_size ->
            reads: [meta, reads]
            k: merylK(genome_size)
        }
    MERYL_COUNT(ch_meryl_in.reads, ch_meryl_in.k)

    // An assembler that produced nothing is not a candidate: every check below would fail
    // on an empty FASTA, and it could not win in any case.
    ch_named = ch_candidates
        .filter { _meta, _assembler, assembly -> assembly.size() > 0 }
        .map { meta, assembler, assembly ->
            [meta + [id: "${meta.id}.${assembler}", sample: meta.id, assembler: assembler], assembly]
        }

    // Trim before measuring: an untrimmed wrap makes an assembly look both longer and
    // less circular than it is, which are two of the things being compared.
    CONTIG_ENDS_CANDIDATE(ch_named)
    CIRCULARISE(ch_named.join(CONTIG_ENDS_CANDIDATE.out.overlaps))

    ch_scored = CIRCULARISE.out.assembly
    ch_candidate_reads = ch_scored.map { meta, _assembly -> [meta.sample, meta] }
        .combine(ch_reads.map { meta, reads -> [meta.id, reads] }, by: 0)
        .map { _sample, meta, reads -> [meta, reads] }

    SEQKIT_STATS_CANDIDATE(ch_scored)
    MERQURY_CANDIDATE(
        ch_scored.map { meta, assembly -> [meta.sample, meta, assembly] }
            .combine(MERYL_COUNT.out.meryl_db.map { meta, db -> [meta.id, db] }, by: 0)
            .map { _sample, meta, assembly, db -> [meta, db, assembly] }
    )
    MAP_READS_CANDIDATE(ch_scored.join(ch_candidate_reads))
    CLIPPING_PILEUPS_CANDIDATE(
        MAP_READS_CANDIDATE.out.bam.join(MAP_READS_CANDIDATE.out.bam_extend)
    )

    // Back to one row per sample: groupTuple waits for every candidate's checks, which is
    // the point of it.
    ch_by_sample = channel.empty()
        .mix(
            SEQKIT_STATS_CANDIDATE.out.stats.map { meta, f -> [meta, 'stats', f] },
            CIRCULARISE.out.circularity.map { meta, f -> [meta, 'circularity', f] },
            MERQURY_CANDIDATE.out.assembly_qv.map { meta, f -> [meta, 'qv', f] },
            MERQURY_CANDIDATE.out.stats.map { meta, f -> [meta, 'completeness', f] },
            MAP_READS_CANDIDATE.out.stats.map { meta, f -> [meta, 'mapping', f] },
            CLIPPING_PILEUPS_CANDIDATE.out.pileups.map { meta, f -> [meta, 'clipping', f] },
        )
        .map { meta, kind, f -> [meta.subMap(meta.keySet() - ['id', 'sample', 'assembler']) + [id: meta.sample], kind, f] }
        .groupTuple()
        .map { meta, kinds, files ->
            def byKind = [stats: [], circularity: [], qv: [], completeness: [], mapping: [], clipping: []]
            [kinds, files].transpose().each { kind, f -> byKind[kind] << f }
            [meta, byKind.stats, byKind.circularity, byKind.qv, byKind.completeness, byKind.mapping, byKind.clipping]
        }

    SCORE_ASSEMBLIES(ch_by_sample.join(ch_genome_size))

    // The candidates as the selector takes them: parallel lists of assembler names and
    // trimmed FASTAs, one pair of lists per sample. `remainder: true` keeps a sample whose
    // candidates all fell over -- it still has the Flye fallback to be finished from.
    ch_candidate_lists = ch_scored
        .map { meta, assembly -> [meta.subMap(meta.keySet() - ['id', 'sample', 'assembler']) + [id: meta.sample], meta.assembler, assembly] }
        .groupTuple()

    SELECT_ASSEMBLY(
        ch_selection_in
            .join(SCORE_ASSEMBLIES.out.scores, remainder: true)
            .join(ch_candidate_lists, remainder: true)
            .filter { row -> row[1] != null }
            .map { meta, consensus, metrics, fallback, scores, assemblers, assemblies ->
                [meta, consensus, metrics, fallback, scores ?: [], assemblers ?: [], assemblies ?: []]
            }
    )

    emit:
    assembly = SELECT_ASSEMBLY.out.assembly
    assembly_source = SELECT_ASSEMBLY.out.summary
    scores = SCORE_ASSEMBLIES.out.scores
    circularity = CIRCULARISE.out.circularity
    meryl_db = MERYL_COUNT.out.meryl_db
}

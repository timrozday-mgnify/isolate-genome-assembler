include { AUTOCYCLER_SUBSAMPLE                        } from '../../modules/local/autocycler_subsample'
include { FLYE                                        } from '../../modules/local/flye'
include { FLYE as FLYE_FULL                           } from '../../modules/local/flye'
include { HIFIASM                                     } from '../../modules/local/hifiasm'
include { HIFIASM as HIFIASM_FULL                     } from '../../modules/local/hifiasm'
include { RAVEN                                       } from '../../modules/local/raven'
include { RAVEN as RAVEN_FULL                         } from '../../modules/local/raven'
include { CANU                                        } from '../../modules/local/canu'
include { CANU as CANU_FULL                           } from '../../modules/local/canu'
include { MINIASM_OVERLAP                             } from '../../modules/local/miniasm_overlap'
include { MINIASM_OVERLAP as MINIASM_OVERLAP_FULL     } from '../../modules/local/miniasm_overlap'
include { MINIASM                                     } from '../../modules/local/miniasm'
include { MINIASM as MINIASM_FULL                     } from '../../modules/local/miniasm'
include { MINIPOLISH                                  } from '../../modules/local/minipolish'
include { MINIPOLISH as MINIPOLISH_FULL               } from '../../modules/local/minipolish'
include { METAMDBG                                    } from '../../modules/local/metamdbg'
include { METAMDBG as METAMDBG_FULL                   } from '../../modules/local/metamdbg'
include { PLASSEMBLER                                 } from '../../modules/local/plassembler'
include { MYLOASM                                     } from '../../modules/local/myloasm'
include { MYLOASM as MYLOASM_FULL                     } from '../../modules/local/myloasm'
include { LJA                                         } from '../../modules/local/lja'
include { LJA as LJA_FULL                             } from '../../modules/local/lja'
include { NORMALISE_HEADERS                           } from '../../modules/local/normalise_headers'
include { NORMALISE_HEADERS as NORMALISE_HEADERS_FULL } from '../../modules/local/normalise_headers'
include { ASSEMBLY_ATTEMPTS                           } from '../../modules/local/assembly_attempts'
include { AUTOCYCLER_CONSENSUS                        } from '../../modules/local/autocycler_consensus'

// Autocycler's subsets are named sample_01.fastq, sample_02.fastq and so on.
def subsetOf(fastq) {
    fastq.simpleName.replaceFirst(/^sample_/, '')
}

// Stages 3 and 4. Subsample, assemble each subset with every enabled assembler, and
// combine the lot into one consensus. Every assembler runs in its own pinned container,
// with the flags Autocycler v0.7.0's helper.rs uses, so each gets the resources it needs
// and one crashing on one subset costs only that assembly.
workflow ASSEMBLY {
    take:
    ch_reads // tuple: [meta, screened HiFi FASTQ]
    ch_read_qc // tuple: [meta, read_qc.tsv], for the genome size the subsampler needs

    main:
    def known = ['flye', 'hifiasm', 'raven', 'canu', 'miniasm', 'metamdbg', 'plassembler', 'myloasm', 'lja']
    def assemblers = params.assemblers.tokenize(',')*.trim()
    def unknown = assemblers - known
    if (unknown) {
        error "Unknown --assemblers value(s): ${unknown.join(', ')}"
    }

    // Plassembler is dropped rather than rejected: it is plasmid-only, so it is never a
    // whole-genome candidate, and it already runs on the full read set in FINISHING.
    def fullAssemblers = (params.full_read_assemblers?.tokenize(',')*.trim() ?: assemblers) - ['plassembler']
    def unknownFull = fullAssemblers - known
    if (unknownFull) {
        error "Unknown --full_read_assemblers value(s): ${unknownFull.join(', ')}"
    }

    ch_genome_size = ch_read_qc
        .splitCsv(header: true, sep: '\t')
        .map { meta, row -> [meta, row.genome_size_used] }

    // A sample with `autocycler_dir` is being reclustered by hand, so it skips the whole
    // expensive half of this subworkflow.
    ch_input = ch_reads
        .join(ch_genome_size)
        .branch { meta, _reads, _genome_size ->
            curated: meta.autocycler_dir
            fresh: true
        }

    AUTOCYCLER_SUBSAMPLE(ch_input.fresh)

    ch_subsets = AUTOCYCLER_SUBSAMPLE.out.subsets
        .transpose()
        .map { meta, genome_size, fastq -> [meta, subsetOf(fastq), fastq, genome_size] }
    ch_assembler_in = ch_subsets.map { meta, subset, fastq, _genome_size -> [meta, subset, fastq] }

    // A disabled assembler is given no work at all, rather than being switched off inside
    // the process, so it leaves no trace of skipped tasks.
    FLYE('flye' in assemblers ? ch_assembler_in : channel.empty())
    HIFIASM('hifiasm' in assemblers ? ch_assembler_in : channel.empty())
    RAVEN('raven' in assemblers ? ch_assembler_in : channel.empty())
    CANU('canu' in assemblers ? ch_subsets : channel.empty())
    METAMDBG('metamdbg' in assemblers ? ch_assembler_in : channel.empty())
    PLASSEMBLER(
        'plassembler' in assemblers ? ch_assembler_in : channel.empty(),
        params.plassembler_db ? file(params.plassembler_db, checkIfExists: true) : [],
    )

    MYLOASM('myloasm' in assemblers ? ch_assembler_in : channel.empty())
    LJA('lja' in assemblers ? ch_assembler_in : channel.empty())

    MINIASM_OVERLAP('miniasm' in assemblers ? ch_assembler_in : channel.empty())
    MINIASM(MINIASM_OVERLAP.out.overlap)
    MINIPOLISH(MINIASM.out.graph)

    ch_native = channel.empty().mix(
        FLYE.out.assembly,
        HIFIASM.out.assembly,
        RAVEN.out.assembly,
        CANU.out.assembly,
        METAMDBG.out.assembly,
        PLASSEMBLER.out.assembly,
        MINIPOLISH.out.assembly,
        MYLOASM.out.assembly,
        LJA.out.assembly,
    )
    NORMALISE_HEADERS(ch_native)

    // Every enabled assembler also runs on the whole read set. These assemblies are scored
    // and one of them is the fallback; they must NOT reach the consensus, because they
    // share reads with every subset and so would cast a vote correlated with all the
    // others.
    ch_full_in = ch_reads.map { meta, reads -> [meta, 'full', reads] }
    ch_full_sized = ch_reads.join(ch_genome_size).map { meta, reads, genome_size -> [meta, 'full', reads, genome_size] }

    // Flye runs unconditionally, whatever --full_read_assemblers says: it is the
    // last-resort fallback, it is cheap next to 28 subset assemblies, and running it
    // always keeps the DAG static and resume-safe.
    FLYE_FULL(ch_full_in)
    HIFIASM_FULL('hifiasm' in fullAssemblers ? ch_full_in : channel.empty())
    RAVEN_FULL('raven' in fullAssemblers ? ch_full_in : channel.empty())
    CANU_FULL('canu' in fullAssemblers ? ch_full_sized : channel.empty())
    METAMDBG_FULL('metamdbg' in fullAssemblers ? ch_full_in : channel.empty())
    MYLOASM_FULL('myloasm' in fullAssemblers ? ch_full_in : channel.empty())
    LJA_FULL('lja' in fullAssemblers ? ch_full_in : channel.empty())

    MINIASM_OVERLAP_FULL('miniasm' in fullAssemblers ? ch_full_in : channel.empty())
    MINIASM_FULL(MINIASM_OVERLAP_FULL.out.overlap)
    MINIPOLISH_FULL(MINIASM_FULL.out.graph)

    ch_native_full = channel.empty().mix(
        FLYE_FULL.out.assembly,
        HIFIASM_FULL.out.assembly,
        RAVEN_FULL.out.assembly,
        CANU_FULL.out.assembly,
        METAMDBG_FULL.out.assembly,
        MINIPOLISH_FULL.out.assembly,
        MYLOASM_FULL.out.assembly,
        LJA_FULL.out.assembly,
    )
    NORMALISE_HEADERS_FULL(ch_native_full)

    // Nothing scores them yet, so selection still takes Flye's.
    ch_full_assemblies = NORMALISE_HEADERS_FULL.out.assembly
        .map { meta, assembly -> [meta, assembly.simpleName.replaceFirst(/_full$/, ''), assembly] }
    ch_flye_full = ch_full_assemblies
        .filter { _meta, assembler, _assembly -> assembler == 'flye' }
        .map { meta, _assembler, assembly -> [meta, assembly] }

    // Every assembly that should have been produced, whether or not it was. Built from the
    // subsets that exist rather than from params, so a subsample that yielded fewer sets
    // than asked for is described honestly.
    ch_expected = ch_subsets
        .flatMap { meta, subset, _fastq, _genome_size -> assemblers.collect { assembler -> [meta, "${assembler}_${subset}"] } }
        .groupTuple()

    // `remainder: true` keeps a sample whose assemblers all failed: it still gets an
    // attempts file, and still reaches the consensus process, which falls back to Flye.
    // An assembler can succeed and still produce nothing -- Plassembler finding no
    // plasmids is the ordinary case -- which leaves an empty FASTA. `autocycler compress`
    // refuses to read an empty input file and aborts, taking the whole consensus with it,
    // so one plasmid-free subset used to cost the sample every other assembly it had.
    // ASSEMBLY_ATTEMPTS still calls the attempt failed: it tests the staged file for
    // content, so a file dropped here and a file never written look the same to it.
    ch_assemblies = NORMALISE_HEADERS.out.assembly
        .filter { _meta, assembly -> assembly.size() > 0 }
        .groupTuple()
    ASSEMBLY_ATTEMPTS(
        ch_expected
            .join(ch_assemblies, remainder: true)
            .map { meta, expected, assemblies -> [meta, expected, assemblies ?: []] }
    )

    ch_consensus_in = ch_input.fresh
        .map { meta, reads, _genome_size -> [meta, reads] }
        .join(ch_assemblies, remainder: true)
        .map { meta, reads, assemblies -> [meta, assemblies ?: [], reads, []] }
        .mix(ch_input.curated.map { meta, reads, _genome_size -> [meta, [], reads, meta.autocycler_dir] })

    AUTOCYCLER_CONSENSUS(ch_consensus_in)

    // Selection happens in SCORING, which needs both the consensus and every scored
    // candidate; what leaves here is the raw material for it.
    ch_selection_in = AUTOCYCLER_CONSENSUS.out.consensus
        .join(AUTOCYCLER_CONSENSUS.out.metrics)
        .join(ch_flye_full, remainder: true)
        // A fallback with no consensus arrives as [meta, null, fallback], one null for
        // the whole missing side, so it must go before the four-way destructuring.
        .filter { row -> row[1] != null }
        .map { meta, consensus, metrics, fallback -> [meta, consensus, metrics, fallback ?: []] }

    emit:
    full_assemblies = ch_full_assemblies
    selection_in = ch_selection_in // tuple: [meta, consensus, consensus metrics, Flye fallback]
    attempts = ASSEMBLY_ATTEMPTS.out.attempts
    autocycler_dir = AUTOCYCLER_CONSENSUS.out.autocycler_dir
    consensus = AUTOCYCLER_CONSENSUS.out.consensus
    consensus_gfa = AUTOCYCLER_CONSENSUS.out.gfa
    autocycler_table = AUTOCYCLER_CONSENSUS.out.table
}

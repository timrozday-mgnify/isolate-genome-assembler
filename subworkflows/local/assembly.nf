include { AUTOCYCLER_SUBSAMPLE                        } from '../../modules/local/autocycler_subsample'
include { FLYE                                        } from '../../modules/local/flye'
include { FLYE as FLYE_FULL                           } from '../../modules/local/flye'
include { HIFIASM                                     } from '../../modules/local/hifiasm'
include { RAVEN                                       } from '../../modules/local/raven'
include { CANU                                        } from '../../modules/local/canu'
include { MINIASM_OVERLAP                             } from '../../modules/local/miniasm_overlap'
include { MINIASM                                     } from '../../modules/local/miniasm'
include { MINIPOLISH                                  } from '../../modules/local/minipolish'
include { METAMDBG                                    } from '../../modules/local/metamdbg'
include { PLASSEMBLER                                 } from '../../modules/local/plassembler'
include { NORMALISE_HEADERS                           } from '../../modules/local/normalise_headers'
include { NORMALISE_HEADERS as NORMALISE_HEADERS_FULL } from '../../modules/local/normalise_headers'
include { ASSEMBLY_ATTEMPTS                           } from '../../modules/local/assembly_attempts'
include { AUTOCYCLER_CONSENSUS                        } from '../../modules/local/autocycler_consensus'
include { SELECT_ASSEMBLY                             } from '../../modules/local/select_assembly'

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
    def assemblers = params.assemblers.tokenize(',')*.trim()
    def unknown = assemblers - ['flye', 'hifiasm', 'raven', 'canu', 'miniasm', 'metamdbg', 'plassembler', 'myloasm', 'lja']
    if (unknown) {
        error "Unknown --assemblers value(s): ${unknown.join(', ')}"
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
    )
    NORMALISE_HEADERS(ch_native)

    // The fallback runs for every sample, unconditionally: it is cheap next to 28 subset
    // assemblies, and running it always keeps the DAG static and resume-safe.
    FLYE_FULL(ch_reads.map { meta, reads -> [meta, 'full', reads] })
    NORMALISE_HEADERS_FULL(FLYE_FULL.out.assembly.map { meta, _assembler, _subset, files -> [meta, 'flye', 'full', files] })

    // Every assembly that should have been produced, whether or not it was. Built from the
    // subsets that exist rather than from params, so a subsample that yielded fewer sets
    // than asked for is described honestly.
    ch_expected = ch_subsets
        .flatMap { meta, subset, _fastq, _genome_size -> assemblers.collect { assembler -> [meta, "${assembler}_${subset}"] } }
        .groupTuple()

    // `remainder: true` keeps a sample whose assemblers all failed: it still gets an
    // attempts file, and still reaches the consensus process, which falls back to Flye.
    ch_assemblies = NORMALISE_HEADERS.out.assembly.groupTuple()
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

    SELECT_ASSEMBLY(
        AUTOCYCLER_CONSENSUS.out.consensus
            .join(AUTOCYCLER_CONSENSUS.out.metrics)
            .join(NORMALISE_HEADERS_FULL.out.assembly, remainder: true)
            .map { meta, consensus, metrics, fallback -> [meta, consensus, metrics, fallback ?: []] }
            .filter { _meta, consensus, _metrics, _fallback -> consensus }
    )

    emit:
    assembly = SELECT_ASSEMBLY.out.assembly
    assembly_source = SELECT_ASSEMBLY.out.summary
    attempts = ASSEMBLY_ATTEMPTS.out.attempts
    autocycler_dir = AUTOCYCLER_CONSENSUS.out.autocycler_dir
    autocycler_table = AUTOCYCLER_CONSENSUS.out.table
}

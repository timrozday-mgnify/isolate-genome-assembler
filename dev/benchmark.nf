// Phase 6 assembler benchmark. Four steps, run from a directory outside the repo:
//
//   1. nextflow run <repo>/dev/benchmark.nf --mode simulate
//        fetches the truth genomes in benchmark_genomes.tsv, simulates HiFi reads for the
//        simulated ones (PBSIM3 CCS mode at --depths), downloads the real runs, and writes
//        benchmark/samples.yml.
//   2. the pipeline itself on benchmark/samples.yml, with every assembler and 6 subsets,
//      publishing the input assemblies:
//        nextflow run <repo>/main.nf --input benchmark/samples.yml --outdir results \
//          --assemblers flye,hifiasm,raven,canu,miniasm,metamdbg,plassembler,myloasm,lja \
//          --subsample_count 6 --publish_input_assemblies <database params>
//   3. nextflow run <repo>/dev/benchmark.nf --mode consensus --results <pipeline outdir>
//        reruns only Autocycler's consensus for each arm in benchmark_arms.tsv, on that
//        arm's subset of the published input assemblies.
//   4. uv run <repo>/dev/assembler_benchmark.py --benchmark benchmark --results <pipeline outdir>
//        scores each arm against the truth genomes and writes dev/assembler_benchmark.{csv,md}.
//
// Arms share one set of input assemblies instead of rerunning the assemblers per arm, so
// the arms differ only in what they choose to combine. Subsets 01-04 of a 6-subset run
// stand in for a 4-subset run: Autocycler sizes each subset from the total depth alone,
// so only which reads overlap between subsets changes.

include { DOWNLOAD_DATABASE as FETCH_GENOME } from '../modules/local/download_database'
include { DOWNLOAD_DATABASE as FETCH_READS  } from '../modules/local/download_database'
include { AUTOCYCLER_CONSENSUS              } from '../modules/local/autocycler_consensus'

// The truth genome as one FASTA, straight from the NCBI Datasets zip.
process TRUTH {
    tag "${name}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/python:3.12.12'
        : 'quay.io/biocontainers/python:3.12.12'}"

    input:
    tuple val(name), path(zip)

    output:
    tuple val(name), path("${name}.fasta"), emit: truth

    script:
    """
    python3 -m zipfile -e ${zip} unzipped
    cat unzipped/ncbi_dataset/data/*/*.fna > ${name}.fasta
    """

    stub:
    """
    printf '>chrom chromosome\\nACGT\\n>plas plasmid p1\\nACGT\\n' > ${name}.fasta
    """
}

// PBSIM3 multi-pass subreads, one replicon at a time. Every replicon is circular, so each is
// simulated from two copies end to end at half the depth, which gives reads across the
// origin; --length-max keeps a read no longer than the replicon itself, up to PBSIM3's own
// 1 Mb limit. Copy number is not modelled: every replicon gets the sample depth.
process SIMULATE_SUBREADS {
    tag "${name}_${depth}x"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/pbsim3:3.0.5--h9948957_2'
        : 'quay.io/biocontainers/pbsim3:3.0.5--h9948957_2'}"

    input:
    tuple val(name), path(truth), val(depth)

    output:
    tuple val(name), val(depth), path('r*_0001.bam'), emit: subreads

    script:
    """
    samtools faidx ${truth}
    i=0
    while read -r contig length _; do
        i=\$((i + 1))
        samtools faidx ${truth} "\$contig" \\
            | awk 'NR == 1 { print ">replicon"; next } { s = s \$0 } END { print s s }' > doubled.fa
        # A replicon shorter than twice the mean read length (small plasmids) gets its mean and
        # sd scaled down to half its length: PBSIM3 dies with SIGFPE when --length-max sits
        # far below --length-mean.
        mean=\$(( length / 2 < ${params.sim_length_mean} ? length / 2 : ${params.sim_length_mean} ))
        sd=\$(( mean * ${params.sim_length_sd} / ${params.sim_length_mean} ))
        pbsim --strategy wgs --method qshmm --qshmm /usr/local/data/QSHMM-RSII.model \\
            --genome doubled.fa --depth ${depth / 2} --pass-num ${params.sim_passes} \\
            --length-mean \$mean --length-sd \$sd \\
            --length-max \$(( length < 1000000 ? length : 1000000 )) --prefix "r\$i" --id-prefix "R\${i}_" --seed \$((${depth} * 100 + i))
        rm doubled.fa "r\${i}_0001.maf.gz" "r\${i}_0001.ref"
    done < ${truth}.fai
    """

    stub:
    """
    touch r1_0001.bam
    """
}

// Subreads to HiFi reads. ccs keeps reads at rq >= 0.99 by default, which is what a
// hifi_reads file holds.
process CCS {
    tag "${name}_${depth}x"
    label 'process_high'

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/pbccs:6.4.0--h9ee0642_0'
        : 'quay.io/biocontainers/pbccs:6.4.0--h9ee0642_0'}"

    input:
    tuple val(name), val(depth), path(subreads)

    output:
    tuple val("${name}_${depth}x"), path("${name}_${depth}x.fastq.gz"), emit: reads

    script:
    """
    for bam in ${subreads}; do
        ccs "\$bam" "\${bam%.bam}.fastq.gz" -j ${task.cpus} --log-level WARN
    done
    cat r*_0001.fastq.gz > ${name}_${depth}x.fastq.gz
    """

    stub:
    """
    echo | gzip -c > ${name}_${depth}x.fastq.gz
    """
}

def tsvRows(path) {
    file(path, checkIfExists: true).splitCsv(header: true, sep: '\t')
}

workflow SIMULATE {
    main:
    def genomes = tsvRows(params.genomes)
    def depths = params.depths.toString().tokenize(',')*.trim()*.toInteger()

    FETCH_GENOME(
        channel.fromList(genomes).map { g ->
            [g.name, "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${g.assembly}/download?include_annotation_type=GENOME_FASTA", "${g.name}.zip"]
        }
    )
    TRUTH(FETCH_GENOME.out.database)

    def simulated = genomes.findAll { g -> !g.reads_url }*.name
    SIMULATE_SUBREADS(TRUTH.out.truth.filter { name, _truth -> name in simulated }.combine(depths))
    CCS(SIMULATE_SUBREADS.out.subreads)

    FETCH_READS(
        channel.fromList(genomes.findAll { g -> g.reads_url }).map { g -> [g.name, g.reads_url, "${g.name}.fastq.gz"] }
    )

    // Paths are relative to samples.yml, which sits next to reads/ and truth/. The truth
    // genome is also the sample's `reference`, so the pipeline's skani and dotplot checks
    // run against it.
    CCS.out.reads
        .map { id, _reads -> [id, id.replaceFirst(/_\d+x$/, '')] }
        .mix(FETCH_READS.out.database.map { name, _reads -> [name, name] })
        .map { id, name -> "- id: ${id}\n  reads: reads/${id}.fastq.gz\n  reference: truth/${name}.fasta\n" }
        .collectFile(name: 'samples.yml', storeDir: params.outdir, sort: true)
}

workflow CONSENSUS {
    main:
    def benchmarkDir = file(params.outdir)
    def samples = new org.yaml.snakeyaml.Yaml().load(benchmarkDir.resolve('samples.yml').text)
    def arms = tsvRows(params.arms)

    ch_consensus_in = channel.fromList(samples).flatMap { sample ->
        def inputs = files("${params.results}/assemblies/${sample.id}/inputs/*.fasta")
        if (!inputs) {
            error "No input assemblies for ${sample.id} under ${params.results}. Run the pipeline with --publish_outputs --publish_input_assemblies first."
        }
        arms.collect { arm ->
            def subsets = (1..arm.subsets.toInteger()).collect { String.format('%02d', it) }
            def wanted = arm.assemblers.tokenize(',').collectMany { a -> subsets.collect { s -> "${a}_${s}.fasta".toString() } }
            [
                [id: "${sample.id}.${arm.arm}", sample: sample.id, arm: arm.arm],
                inputs.findAll { it.name in wanted },
                benchmarkDir.resolve(sample.reads),
                [],
            ]
        }
    }
    AUTOCYCLER_CONSENSUS(ch_consensus_in)
}

workflow {
    if (params.mode == 'simulate') {
        SIMULATE()
    }
    else if (params.mode == 'consensus') {
        CONSENSUS()
    }
    else {
        error "--mode must be 'simulate' or 'consensus', not '${params.mode}'"
    }
}

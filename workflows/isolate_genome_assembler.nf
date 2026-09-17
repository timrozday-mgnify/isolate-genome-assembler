workflow ISOLATE_GENOME_ASSEMBLER {
    take:
    ch_samples // tuple: [meta map, one or more PacBio HiFi read files]

    main:
    ch_samples.view { meta, reads -> "Validated sample ${meta.id}: ${reads.size()} read input(s)" }

    emit:
    samples = ch_samples
}

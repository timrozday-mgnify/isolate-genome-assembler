include { DOWNLOAD_DATABASE     } from '../../modules/local/download_database'
include { SYLPH_SKETCH_GENOMES  } from '../../modules/local/sylph_sketch_genomes'
include { PLASSEMBLER_DOWNLOAD  } from '../../modules/local/plassembler_download'

// `--prepare_databases`. Run once per site into a shared location; the pipeline
// proper only reads the resulting paths through the database params.
workflow PREPARE_DATABASES {
    main:
    ch_downloads = Channel.fromList([
        ['sylph_gtdb_db', params.sylph_gtdb_db_url, file(params.sylph_gtdb_db_url).name],
        ['sylph_gtdb_taxonomy', params.sylph_gtdb_taxonomy_url, file(params.sylph_gtdb_taxonomy_url).name],
        ['human_chm13', params.chm13_url, file(params.chm13_url).name],
        ['human_grch38', params.grch38_url, file(params.grch38_url).name],
    ])

    DOWNLOAD_DATABASE(ch_downloads)

    ch_human_genomes = DOWNLOAD_DATABASE.out.database
        .filter { name, _genome -> name in ['human_chm13', 'human_grch38'] }
        .map { _name, genome -> genome }
        .collect()
        .map { genomes -> ['sylph_human_db', genomes] }

    SYLPH_SKETCH_GENOMES(ch_human_genomes)
    PLASSEMBLER_DOWNLOAD()

    emit:
    databases = DOWNLOAD_DATABASE.out.database
        .mix(SYLPH_SKETCH_GENOMES.out.database)
        .mix(PLASSEMBLER_DOWNLOAD.out.database)
}

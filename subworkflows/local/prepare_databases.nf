include { DOWNLOAD_DATABASE     } from '../../modules/local/download_database'
include { SYLPH_SKETCH_GENOMES  } from '../../modules/local/sylph_sketch_genomes'
include { PLASSEMBLER_DOWNLOAD  } from '../../modules/local/plassembler_download'
include { UNTAR                 } from '../../modules/nf-core/untar'
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload'
include { BUSCO_DOWNLOAD        } from '../../modules/nf-core/busco/download'
include { DIAMOND_MAKEDB        } from '../../modules/nf-core/diamond/makedb'

// `--prepare_databases`. Run once per site into a shared location; the pipeline
// proper only reads the resulting paths through the database params. A database whose
// param is already set is skipped, so a partial build can be finished without refetching.
def missing(name) {
    !(params[name] != null && params[name].toString().trim())
}

workflow PREPARE_DATABASES {
    main:
    // CHM13 is both --human_reference and half of the human sylph sketch.
    def downloads = [
        ['sylph_gtdb_db', params.sylph_gtdb_db_url, missing('sylph_gtdb_db')],
        ['sylph_gtdb_taxonomy', params.sylph_gtdb_taxonomy_url, missing('sylph_gtdb_taxonomy')],
        ['human_chm13', params.chm13_url, missing('human_reference') || missing('sylph_human_db')],
        ['human_grch38', params.grch38_url, missing('sylph_human_db')],
        ['checkm2_db', params.checkm2_db_url, missing('checkm2_db')],
        ['gtdbtk_db', params.gtdbtk_db_url, missing('gtdbtk_db')],
        ['ideel_db', params.ideel_db_url, missing('ideel_db')],
    ]
        .findAll { _name, _url, wanted -> wanted }
        .collect { name, url, _wanted -> [name, url, file(url).name] }

    DOWNLOAD_DATABASE(Channel.fromList(downloads))
    ch_databases = DOWNLOAD_DATABASE.out.database

    if (missing('sylph_human_db')) {
        ch_human_genomes = DOWNLOAD_DATABASE.out.database
            .filter { name, _genome -> name in ['human_chm13', 'human_grch38'] }
            .map { _name, genome -> genome }
            .collect()
            .map { genomes -> ['sylph_human_db', genomes] }
        SYLPH_SKETCH_GENOMES(ch_human_genomes)
        ch_databases = ch_databases.mix(SYLPH_SKETCH_GENOMES.out.database)
    }

    if (missing('plassembler_db')) {
        PLASSEMBLER_DOWNLOAD()
        ch_databases = ch_databases.mix(PLASSEMBLER_DOWNLOAD.out.database)
    }

    UNTAR(
        DOWNLOAD_DATABASE.out.database
            .filter { name, _archive -> name in ['checkm2_db', 'gtdbtk_db'] }
            .map { name, archive -> [[id: name], archive] }
    )
    DIAMOND_MAKEDB(
        DOWNLOAD_DATABASE.out.database
            .filter { name, _fasta -> name == 'ideel_db' }
            .map { name, fasta -> [[id: name], fasta] },
        [],
        [],
        [],
    )
    ch_databases = ch_databases
        .mix(UNTAR.out.untar.map { meta, dir -> [meta.id, dir] })
        .mix(DIAMOND_MAKEDB.out.db.map { meta, db -> [meta.id, db] })

    if (missing('bakta_db')) {
        BAKTA_BAKTADBDOWNLOAD()
        ch_databases = ch_databases.mix(BAKTA_BAKTADBDOWNLOAD.out.db.map { db -> ['bakta_db', db] })
    }

    if (missing('busco_db')) {
        BUSCO_DOWNLOAD(params.busco_lineage)
        ch_databases = ch_databases.mix(BUSCO_DOWNLOAD.out.download_dir.map { dir -> ['busco_db', dir] })
    }

    emit:
    databases = ch_databases
}

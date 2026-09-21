include { ISOLATE_GENOME_ASSEMBLER } from './workflows/isolate_genome_assembler'
include { PREPARE_DATABASES as BUILD_DATABASES } from './subworkflows/local/prepare_databases'

def isNonBlank(value) {
    value != null && value.toString().trim()
}

def failSamplesheet(message) {
    error "Invalid samplesheet '${params.input}': ${message}"
}

def resolveInputFile(value, baseDir, sampleId, field, extensions) {
    if (!(value instanceof String) || !isNonBlank(value)) {
        failSamplesheet("sample '${sampleId}' key '${field}' must be a non-empty path")
    }

    def candidate = java.nio.file.Path.of(value).normalize()
    def resolved = candidate.isAbsolute() ? candidate : baseDir.resolve(candidate).normalize()
    if (!java.nio.file.Files.isRegularFile(resolved)) {
        failSamplesheet("sample '${sampleId}' key '${field}' does not exist or is not a file: ${resolved}")
    }
    if (!(resolved.fileName.toString() ==~ extensions)) {
        failSamplesheet("sample '${sampleId}' key '${field}' has an unsupported file type: ${resolved.fileName}")
    }
    file(resolved)
}

def resolveInputDirectory(value, baseDir, sampleId, field) {
    if (!(value instanceof String) || !isNonBlank(value)) {
        failSamplesheet("sample '${sampleId}' key '${field}' must be a non-empty directory path")
    }

    def candidate = java.nio.file.Path.of(value).normalize()
    def resolved = candidate.isAbsolute() ? candidate : baseDir.resolve(candidate).normalize()
    if (!java.nio.file.Files.isDirectory(resolved)) {
        failSamplesheet("sample '${sampleId}' key '${field}' does not exist or is not a directory: ${resolved}")
    }
    file(resolved)
}

def parseSamplesheet(samplesheet) {
    def samplesheetPath = file(samplesheet, checkIfExists: true)
    def document
    try {
        document = new org.yaml.snakeyaml.Yaml().load(samplesheetPath.text)
    } catch (Exception exception) {
        failSamplesheet("could not parse YAML: ${exception.message}")
    }

    def rows
    if (document instanceof List) {
        rows = document
    } else if (document instanceof Map && document.keySet() == ['samples'] as Set) {
        if (document.samples instanceof List) {
            rows = document.samples
        } else if (document.samples instanceof Map) {
            rows = document.samples.collect { sampleId, values ->
                if (!(values instanceof Map)) {
                    failSamplesheet("sample '${sampleId}' must be a mapping")
                }
                if (values.id != null && values.id.toString() != sampleId.toString()) {
                    failSamplesheet("sample map key '${sampleId}' conflicts with id '${values.id}'")
                }
                values + [id: sampleId.toString()]
            }
        } else {
            failSamplesheet("the 'samples' value must be a list or an id-keyed mapping")
        }
    } else {
        failSamplesheet("top level must be a list or a mapping with only a 'samples' key")
    }

    def allowedKeys = [
        'id',
        'reads',
        'genome_size',
        'expected_taxon',
        'reference',
        'autocycler_cluster_args',
        'autocycler_dir'
    ] as Set
    def sampleIds = [] as Set

    rows.withIndex().collect { row, index ->
        if (!(row instanceof Map)) {
            failSamplesheet("entry ${index + 1} must be a mapping")
        }
        def unknownKeys = row.keySet().collect { it.toString() }.findAll { !(it in allowedKeys) }
        if (unknownKeys) {
            failSamplesheet("sample entry ${index + 1} has unsupported key(s): ${unknownKeys.join(', ')}")
        }
        if (!(row.id instanceof String) || !(row.id ==~ /^[A-Za-z0-9][A-Za-z0-9_.-]*$/)) {
            failSamplesheet("entry ${index + 1} key 'id' must contain only letters, digits, '.', '_' or '-'")
        }
        if (!sampleIds.add(row.id)) {
            failSamplesheet("duplicate sample id '${row.id}'")
        }
        if (row.reads == null) {
            failSamplesheet("sample '${row.id}' is missing required key 'reads'")
        }

        def readValues = row.reads instanceof List ? row.reads : [row.reads]
        if (readValues.isEmpty()) {
            failSamplesheet("sample '${row.id}' key 'reads' must contain at least one path")
        }
        def readFiles = readValues.collect { read ->
            resolveInputFile(read, samplesheetPath.parent, row.id, 'reads', ~/(?i).*\.(fastq|fq)(\.gz)?|.*\.bam/)
        }

        def meta = [id: row.id]
        if (isNonBlank(row.genome_size)) {
            if (!(row.genome_size.toString() ==~ /(?i)^\d+(\.\d+)?[kmg]?$/)) {
                failSamplesheet("sample '${row.id}' key 'genome_size' must be a positive size such as '5.2m'")
            }
            meta.genome_size = row.genome_size.toString()
        }
        if (isNonBlank(row.expected_taxon)) {
            meta.expected_taxon = row.expected_taxon.toString()
        }
        if (row.reference != null) {
            meta.reference = resolveInputFile(row.reference, samplesheetPath.parent, row.id, 'reference', ~/(?i).*\.f(ast)?a(\.gz)?/)
        }
        if (isNonBlank(row.autocycler_cluster_args)) {
            meta.autocycler_cluster_args = row.autocycler_cluster_args.toString()
        }
        if (row.autocycler_dir != null) {
            meta.autocycler_dir = resolveInputDirectory(row.autocycler_dir, samplesheetPath.parent, row.id, 'autocycler_dir')
        }

        [meta, readFiles]
    }
}

def requireDatabases() {
    def missing = [
        'sylph_gtdb_db',
        'sylph_gtdb_taxonomy',
        'sylph_human_db',
        'human_reference',
        'plassembler_db',
        'checkm2_db',
        'bakta_db',
        'busco_db',
        'gtdbtk_db',
        'ideel_db',
    ]
        .findAll { !isNonBlank(params[it]) }
    if (missing) {
        error "Missing required database param(s): ${missing.collect { "--${it}" }.join(', ')}. Build them once with `nextflow run main.nf --prepare_databases`."
    }
}

workflow {
    // Run once per site: `nextflow run main.nf --prepare_databases --database_dir /shared/dbs`.
    // Nextflow's strict syntax dropped `-entry`, so the database build is a param, not an
    // entry workflow.
    if (params.prepare_databases) {
        BUILD_DATABASES()
    }
    else {
        ch_samples = params.input ? Channel.fromList(parseSamplesheet(params.input)) : Channel.empty()

        if (params.input) {
            requireDatabases()
        }
        else {
            log.info 'No --input supplied: validated an empty workflow. See assets/samplesheet.example.yml for the required YAML format.'
        }

        ISOLATE_GENOME_ASSEMBLER(ch_samples)
    }
}

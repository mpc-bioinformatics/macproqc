include { GUNZIP } from '../../../modules/nf-core/gunzip/main'
include { TDF2MZML } from '../../../modules/local/tdf2mzml/main'
include { THERMORAWFILEPARSER } from '../../../modules/nf-core/thermorawfileparser/main'
include { UNTAR } from '../../../modules/nf-core/untar/main'
include { UNZIP } from '../../../modules/nf-core/unzip/main'

workflow PREPARE_SPECTRA {
    take:
    ch_spectra // channel: [ val(meta), [ .d-folder | .raw | .mzml ] ]

    main:

    def VENDOR_BRUKER = 'bruker'
    def VENDOR_THERMO = 'thermo'
    def VENDOR_MZML = 'mzml'
    def VENDOR_UNKNOWN = 'unknown'

    def COMPRESSION_GZ = 'gz'
    def COMPRESSION_TAR = 'tar'
    def COMPRESSION_ZIP = 'zip'
    def COMPRESSION_NONE = 'none'
    def COMPRESSION_UNSUPPORTED = 'unsupported'

    // map to assign vendor based on the file extension
    def EXTENSION_TO_VENDOR_MAP = [
        '.d': VENDOR_BRUKER,
        '.d.tar.gz': VENDOR_BRUKER,
        '.d.tar.bz2': VENDOR_BRUKER,
        '.d.zip': VENDOR_BRUKER,
        '.raw': VENDOR_THERMO,
        '.raw.zip': VENDOR_THERMO,
        '.raw.tar.gz': VENDOR_THERMO,
        '.raw.tar.bz2': VENDOR_THERMO,
        '.raw.gz': VENDOR_THERMO,
        '.mzml': VENDOR_MZML,
        '.mzml.zip': VENDOR_MZML,
        '.mzml.tar.gz': VENDOR_MZML,
        '.mzml.tar.bz2': VENDOR_MZML,
        '.mzml.gz': VENDOR_MZML,
    ].asImmutable()

    // map to assign compression type based on the file extension
    def EXTENSION_TO_COMPRESSION_MAP = [
        '.d': COMPRESSION_NONE,
        '.d.tar.gz': COMPRESSION_TAR,
        '.d.tar.bz2': COMPRESSION_TAR,
        '.d.zip': COMPRESSION_ZIP,
        '.raw': COMPRESSION_NONE,
        '.raw.zip': COMPRESSION_ZIP,
        '.raw.tar.gz': COMPRESSION_TAR,
        '.raw.tar.bz2': COMPRESSION_TAR,
        '.raw.gz': COMPRESSION_GZ,
        '.mzml': COMPRESSION_NONE,
        '.mzml.zip': COMPRESSION_ZIP,
        '.mzml.tar.gz': COMPRESSION_TAR,
        '.mzml.tar.bz2': COMPRESSION_TAR,
        '.mzml.gz': COMPRESSION_GZ,
    ].asImmutable()

    // add the vendor information to the metadata based on the file extension, if the extension is unknown assign VENDOR_UNKNOWN
    ch_spectra = ch_spectra.map { meta, path ->
        def name = path.name.toLowerCase()
        // longest known extension the file name ends with (file names may contain further dots)
        def full_extension = EXTENSION_TO_VENDOR_MAP.keySet()
            .findAll { ext -> name.endsWith(ext) }
            .max { ext -> ext.length() }
        meta.vendor = EXTENSION_TO_VENDOR_MAP[full_extension] ?: VENDOR_UNKNOWN
        meta.compression = EXTENSION_TO_COMPRESSION_MAP[full_extension] ?: COMPRESSION_UNSUPPORTED
        return [meta, path]
    }

    // Filter the spectra files which are compresseds
    ch_spectra
        .branch { item ->
            gzipped: item[0].compression == COMPRESSION_GZ
            tarred: item[0].compression == COMPRESSION_TAR
            zipped: item[0].compression == COMPRESSION_ZIP
            uncompressed: item[0].compression == COMPRESSION_NONE
            unsupported: true
        }
        .set { ch_branched_spectra }

    ch_branched_spectra.unsupported.map { item ->
        log.warn("Found spectra files with unsupported format: `${item[1]}`. These will be ignored.")
    }

    ch_uncompressed_spectra = ch_branched_spectra.uncompressed

    // Uncompress the spectra files
    GUNZIP(ch_branched_spectra.gzipped)
    ch_gunzipped_spectra = GUNZIP.out.gunzip.map { meta, path -> path.isFile() ? [meta, path] : [meta, path.listFiles().first()] }
    ch_uncompressed_spectra = ch_uncompressed_spectra.mix(ch_gunzipped_spectra)

    UNTAR(ch_branched_spectra.tarred)
    ch_untarred_spectra = UNTAR.out.untar.map { meta, path -> path.isFile() || meta.vendor == VENDOR_BRUKER ? [meta, path] : [meta, path.listFiles().first()] }
    ch_uncompressed_spectra = ch_uncompressed_spectra.mix(ch_untarred_spectra)

    UNZIP(ch_branched_spectra.zipped)
    unzipped_spectra = UNZIP.out.unzipped_archive.map { meta, path -> path.isFile() ? [meta, path] : [meta, path.listFiles().first()] }
    ch_uncompressed_spectra = ch_uncompressed_spectra.mix(unzipped_spectra)

    // remove compression from meta
    ch_uncompressed_spectra = ch_uncompressed_spectra.map { meta, path ->
        meta.remove('compression')
        return [meta, path]
    }

    // Branch into the different kind of spectra files and convert to mzML if necessary
    ch_uncompressed_spectra
        .branch { item ->
            brukerd: item[0].vendor == VENDOR_BRUKER
            thermoraw: item[0].vendor == VENDOR_THERMO
            mzml: item[0].vendor == VENDOR_MZML
        }
        .set { ch_branched_uncompressed_spectra }

    ch_mzmls = ch_branched_uncompressed_spectra.mzml

    TDF2MZML(ch_branched_uncompressed_spectra.brukerd)
    ch_mzmls = ch_mzmls.mix(TDF2MZML.out.spectra)

    THERMORAWFILEPARSER(ch_branched_uncompressed_spectra.thermoraw)
    ch_mzmls = ch_mzmls.mix(THERMORAWFILEPARSER.out.spectra)

    emit:
    mzmls = ch_mzmls
    uncompressed = ch_uncompressed_spectra
}

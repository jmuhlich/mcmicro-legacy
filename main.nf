#!/usr/bin/env nextflow


include {staging}        from './modules/staging'
include {illumination}   from './modules/illumination'
include {registration}   from './modules/registration'
include {dearray}        from './modules/dearray'
include {segmentation}   from './modules/segmentation'
include {quantification} from './modules/quantification'
include {downstream}     from './modules/downstream'
include {viz}            from './modules/viz'
include {background}     from './modules/background'

// Helper functions for finding raw images and precomputed intermediates
def findFiles0(pre, key, pattern) {
    pre[key] ? Channel.fromPath("${params.in}/$key/$pattern") : Channel.empty()
}
def findFiles(pre, key, pattern, ife) {
    pre[key] ? Channel.fromPath("${params.in}/$key/$pattern").ifEmpty(ife) : Channel.empty()
}
def findDirs(pre, key, ife) {
    pre[key] ? Channel.fromPath("${params.in}/$key/*", type: 'dir').ifEmpty(ife) : Channel.empty()
}

// Define the primary mcmicro workflow
workflow {

    main:

    // Expecting --in parameter
    if( !params.containsKey('in') )
        error "Please specify the project directory with --in"

    // Parse MCMICRO parameters (mcp)
    def mcp = mcmicro.Opts.parseParams(
        params, 
        "$projectDir/config/schema.yml",
        "$projectDir/config/defaults.yml"
    )

    // Separate out workflow parameters (wfp) to simplify code
    def wfp = mcp.workflow

    // Identify relevant precomputed intermediates
    // The actual paths to intermediate files are given by
    //   pre.collect{ "${params.in}/$it" }
    def pre = mcmicro.Flow.precomputed(wfp)

    // Check that deprecated locations are empty
    Channel.fromPath( "${params.in}/illumination_profiles/*" )
        .subscribe{ it ->
        error "illumination_profiles/ is deprecated; please use illumination/ instead"
    }

    // Identify marker information
    def chMrk = Channel.fromPath( "${params.in}/markers.csv", checkIfExists: true )

    // Some image formats store multiple fields of view in a single file. Other
    // formats store each field separately, typically in .tif files, with a separate
    // index file to tie them together. We will look for the index files from
    // multiple-file formats in a first, separate pass in order to avoid finding the
    // individual .tif files instead. If no multi-file formats are detected, then we
    // look for the single-file formats. Also, for multi-file formats we need to
    // stage the parent directory and not just the index file.
    def (formatType, formatPattern) =
        files("${params.in}/raw/**${wfp['multi-formats']}") ?
        ["multi", wfp['multi-formats']] : ["single", wfp['single-formats']]

    def stagingDirs = findDirs(pre, 'staging', 
        {error "No subdirectories found in staging directory"})
    def staging_in = stagingDirs
        .map{ tuple(
            mcmicro.Util.getSampleName(it, file("${params.in}/staging")),
            mcmicro.Util.getCycleNameFromDir(it, file("${params.in}/staging")),
            formatType == "single" ? it : it.parent
        )}
    // Here we assemble tuples of 1) path to stage for each raw image (might be a
    // directory) and 2) relative path to the main file for each image. Processes
    // must input the first as a path and the second as a val to avoid incorrect or
    // redundant file staging. They must also only use the second (relative) path to
    // construct pathnames for scripts etc. mcmicro.Util.escapePathForShell must be
    // used when interpolating these paths into script strings, as we are bypassing
    // the normal way that paths are passed to channels which handles this escaping
    // automatically.
    def rawFiles = findFiles(pre, 'raw', "**${formatPattern}",
                         {error "No images found in ${params.in}/raw"})
    def raw = rawFiles
        .map{ tuple(
            mcmicro.Util.getSampleName(it, file("${params.in}/raw")),
            formatType == "single" ? it : it.parent, 
            it
        )}
        .map{ sampleName, toStage, relPath -> 
            tuple(sampleName, toStage, toStage.parent.relativize(relPath).toString()) }

    // Find precomputed intermediates
    def pre_dfp   = findFiles0(pre, 'illumination', "**-dfp.tif")
        .map{ tuple(mcmicro.Util.getSampleName(it, file("${params.in}/illumination")), it) }
    def pre_ffp   = findFiles0(pre, 'illumination', "**-ffp.tif")
        .map{ tuple(mcmicro.Util.getSampleName(it, file("${params.in}/illumination")), it) }
    def pre_img   = findFiles(pre, 'registration', "*.{ome.tiff,ome.tif,tif,tiff,btf,qptiff}",
        {error "No pre-stitched image in ${params.in}/registration"})
    def pre_bsub  = findFiles(pre, 'background', "*.ome.tif",
        {error "No background subtracted image in ${params.in}/background"})
    def pre_bsubm = findFiles(pre, 'background', "*.csv",
        {error "No background subtracted markers file in ${params.in}/background"})
    def pre_cores = findFiles(pre, 'dearray', "*.tif",
        {error "No TMA cores in ${params.in}/dearray"})
    def pre_masks = findFiles(pre, 'dearray', "masks/*.tif",
        {error "No TMA masks in ${params.in}/dearray/masks"})
    def pre_pmap  = findFiles(pre, 'probability-maps', "*/*-pmap.tif",
        {error "No probability maps found in ${params.in}/probability-maps"})
        .map{ f -> tuple(f.getParent().getName(), f) }
        .filter{ wfp['segmentation'].contains(it[0]) }
    def pre_seg   = findFiles(pre, 'segmentation', "**.tif",
        {error "No segmentation masks in ${params.in}/segmentation"})
        .map{ f -> tuple(f.getParent().getName(), f) }.groupTuple()
    def pre_qty   = findFiles(pre, 'quantification', "*.csv",
        {error "No quantification tables in ${params.in}/quantification"})

    if (mcmicro.Flow.doirun('staging', wfp)) {
        staging(mcp, staging_in, chMrk)
        staging.out.map{
            sample, cycle, path ->
            tuple(sample, cycle, path, path.toString().split('/').last())
        }.toSortedList { a, b -> a[1] <=> b[1] }
            .flatMap()
            .map{
                sample, cycle, path, name ->
                tuple(sample, path, name)
            }.set{ sorted_staging }
        raw = raw.mix(sorted_staging)
    }

    def ffp = pre_ffp
    def dfp = pre_dfp
    if (mcmicro.Flow.doirun('illumination', wfp)) {
        illumination(wfp, mcp.modules['illumination'], raw)
        ffp = illumination.out.ffp.mix(ffp)
        dfp = illumination.out.ffp.mix(dfp)
    }
    def img = pre_img
    if (mcmicro.Flow.doirun('registration', wfp)) {
        registration(mcp, raw, ffp, dfp)
        img = registration.out.mix(img)
    }

    // Should background subtraction be applied?
    img = img.branch{
        nobs: !wfp.background
        bs: wfp.background
    }
    chMrk = chMrk.branch{
        nobs: !wfp.background
        bs: wfp.background
    }
    def bsub_image = pre_bsub
    def bsub_marker = pre_bsubm
    // Apply background if specified
    if (mcp.workflow["background"]) {
        background(mcp, img.bs, chMrk.bs)
        bsub_image = background.out.image.mix(bsub_image)
        bsub_marker = background.out.marker.mix(bsub_marker)
    }
    // Reconcile non-background subtracted and background 
    // subtracted images for downstream processing
    img = img.nobs.mix(bsub_image)
    // Reconcile the marker file to the background subtracted csv
    chMrk = chMrk.nobs.mix(bsub_marker)

    // Are we working with a TMA or a whole-slide image?
    img = img.branch{
        wsi: !wfp.tma
        tma: wfp.tma
    }

    def tmacores = pre_cores
    def tmamasks = pre_masks
    // Apply dearray to TMAs only
    if (mcmicro.Flow.doirun('dearray', wfp)) {
        dearray(mcp, img.tma)
        tmacores = dearray.out.cores.mix(tmacores)
        tmamasks = dearray.out.masks.mix(tmamasks)
    }

    // Reconcile WSI and TMA processing for downstream segmentation
    def allimg = img.wsi.mix(tmacores)
    segmentation(mcp, allimg, tmamasks, pre_pmap)

    // Merge segmentation masks against precomputed ones and append markers.csv
    def segMsk = segmentation.out.mix(pre_seg)
    def sft = pre_qty
    if (mcmicro.Flow.doirun('quantification', wfp)) {
        quantification(mcp, allimg, segMsk, chMrk)
        sft = quantification.out.mix(sft)
    }

    // Spatial feature tables -> cell state calling
    if (mcmicro.Flow.doirun('downstream', wfp)) {
        downstream(mcp, sft)
    }

    // Vizualization
    if (mcmicro.Flow.doirun('viz', wfp)) {
        viz(mcp, allimg, chMrk)
    }

    onComplete:

    // Write out parameters used
    def path_qc = "${params.in}/qc"
    
    // Create a provenance directory
    file(path_qc).mkdirs()

    // Write out MCMICRO parameters
    def style = new org.yaml.snakeyaml.DumperOptions();
    style.setPrettyFlow(true);
    style.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK);
    file("${path_qc}/params.yml").withWriter{ out -> 
        new org.yaml.snakeyaml.Yaml(style).dump(mcp, out) 
    }

    // Store additional metadata
    file("${path_qc}/metadata.yml").withWriter{ out ->
        out.println("githubTag: $workflow.revision")
        out.println("githubCommit: $workflow.commitId")
        out.println("roadie: $params.roadie")
    }

}

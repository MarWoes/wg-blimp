
## HACK HACK HACK HACK HACK
# shiny does not support range byte request, so we hack the static file handling function to support range byte requests.
# see https://github.com/rstudio/shiny/pull/747/files
# also, shiny does not support symlinking, so we manually add that as well.
# see https://github.com/rstudio/shiny/issues/1064
hackedStaticHandling <- function(routeMap) {

  return(function(req) {
    if (!identical(req$REQUEST_METHOD, 'GET'))
      return(NULL)

    if (!req$PATH_INFO %in% names(routeMap))
      return(NULL)

    path <- routeMap[[req$PATH_INFO]]

    # shiny doesn't correctly handle symlinks, so we fix those on our own....
    symlinkPath <- Sys.readlink(path)

    if (file.exists(symlinkPath)) {
      abs.path <- symlinkPath
    } else {
      abs.path <- path
    }

    if (is.null(abs.path))
      return(NULL)

    content.type <- shiny:::getContentType(abs.path)
    content.length <- file.info(abs.path)$size

    if(length(req$HTTP_RANGE) && grepl('^bytes=.+', req$HTTP_RANGE)) {
      rng <- as.numeric(
        strsplit(gsub('^bytes=', '', req$HTTP_RANGE), '-')[[1]]
      )
      file.connection <- file(abs.path, "rb")
      seek(file.connection, where = rng[1], origin = "start")
      response.content <- readBin(
        file.connection, 'raw', n=length(rng[1]:rng[2])
      )
      close(file.connection)
      response.headers <- list(`Accept-Ranges` = 'bytes', `Content-Range` = paste("bytes ", rng[1], "-", rng[2], "/", content.length, sep = ""))
      response.code <- 206
    } else {

      # shiny doesn't handle direct download of big files well, so we just don't allow it
      # to prevent someone accidentally downloading such a file and thus crashing everything
      bytes_in_one_gigabyte <- 10e9
      if (content.length >= bytes_in_one_gigabyte)
        return(shiny:::httpResponse(400, content="<h1>File too big for whole download. Use range requests instead.</h1>"))

      response.content <- readBin(abs.path, 'raw', n=content.length)
      response.headers <- list(`Accept-Ranges` = 'bytes')
      response.code <- 200
    }
    return(shiny:::httpResponse(response.code, content.type, response.content, headers = response.headers))
  })
}

## END HACK HACK HACK HACK HACK

shiny.wgbs.serveBamFiles <- function (datasets) {

  servedFilesSubPath <- lapply(datasets, function (dataset) list.files(dataset[["bamDirectory"]], pattern = "bai$|bam$"))
  servedFilesFullPath <- lapply(datasets, function (dataset) list.files(dataset[["bamDirectory"]], pattern = "bai$|bam$", full.names = TRUE))

  datasetLengths <- sapply(names(servedFilesSubPath), function (datasetName) length(servedFilesSubPath[[datasetName]]))

  routeKeys <- paste0("/", rep(names(servedFilesSubPath), datasetLengths), "/", unlist(servedFilesSubPath))

  routeMap <- setNames(unlist(servedFilesFullPath), routeKeys)

  shiny:::handlerManager$addHandler(hackedStaticHandling(routeMap), key = "alignments")

}


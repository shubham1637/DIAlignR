#' Subset an XIC file
#'
#' Create a sqMass file that has chromatograms for given native IDs.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2022) + GPL-3
#' Date: 2022-04-19
#' @param nativeIDs (integer) transition IDs to be kept.
#' @param xicFileIn (character) name to the current file.
#' @param xicFileOut (character) name of the new file.
#' @return (None)
#' @seealso \code{\link{createSqMass}, \link{getNativeIDs}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' oswIn <- file.path(dataPath, "osw", "merged.osw")
#' params <- paramsDIAlignR()
#' params[["context"]] <- "experiment-wide"
#' ids <- getNativeIDs(oswIn, 1338L, params)
#' xicFileIn <- file.path(dataPath, "xics", "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.sqMass")
#' reduceXICs(ids, xicFileIn, xicFileOut = "temp.chrom.sqMass")
#' file.remove("temp.chrom.sqMass")
#' @export
reduceXICs <- function(nativeIDs, xicFileIn, xicFileOut){
  conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = xicFileIn)
  chromDf <- DBI::dbReadTable(conn, "CHROMATOGRAM")
  dataDf <- DBI::dbReadTable(conn, "DATA")
  DBI::dbDisconnect(conn); rm(conn)
  chromDf <- chromDf[chromDf$NATIVE_ID %in% nativeIDs,]
  oldIDs <- chromDf$ID
  newIDs <- seq(from=0, to=(nrow(chromDf)-1), by = 1L)
  dataDf <- dataDf[dataDf$CHROMATOGRAM_ID %in% chromDf$ID,]
  dataDf$CHROMATOGRAM_ID <- newIDs[match(dataDf$CHROMATOGRAM_ID, oldIDs)]
  chromDf$ID <- newIDs
  rownames(dataDf) <- NULL
  rownames(chromDf) <- NULL
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname=":memory:")
  DBI::dbWriteTable(conn=con, name="DATA", dataDf, append=TRUE, row.names = FALSE)
  DBI::dbWriteTable(conn=con, name="CHROMATOGRAM", chromDf, append=TRUE, row.names = FALSE)
  DBI::dbExecute(con, "CREATE INDEX data_chr_idx ON DATA(CHROMATOGRAM_ID);")
  tryCatch(expr = {db <- DBI::dbConnect(RSQLite::SQLite(), dbname=xicFileOut)
  RSQLite::sqliteCopyDatabase(con, db)},
  finally = DBI::dbDisconnect(db))
  DBI::dbDisconnect(con)
  message("Written file ", xicFileOut)
}


#' Fetch NativeIDs
#'
#' Get transition(native) IDs for the peptides.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2022) + GPL-3
#' Date: 2022-04-19
#'
#' @inheritParams alignTargetedRuns
#' @param oswIn (string) path to the osw feature file.
#' @return (integer) a vector of transition IDs.
#' @seealso \code{\link{reduceXICs}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' oswIn <- file.path(dataPath, "osw", "merged.osw")
#' peps <- c(3L, 11L) # No transitions for 3L.
#' params <- paramsDIAlignR()
#' params[["context"]] <- "experiment-wide"
#' getNativeIDs(oswIn, peps, params) # 106468 106469 106470 106471 106472 106473
#' @export
getNativeIDs <- function(oswIn, peps, params=paramsDIAlignR(), oswMerged=TRUE){
  fileInfo <- data.frame("featureFile" = oswIn)
  start_time <- Sys.time()
  precursors <- getPrecursors(fileInfo, oswMerged, params[["runType"]], params[["context"]], params[["maxPeptideFdr"]], params[["level"]])
  if(!is.null(peps)){
    precursors <- precursors[peptide_id %in% peps, ]
    if(nrow(precursors) == 0L) stop("No peptide IDs are found in osw files.")
    data.table::setkeyv(precursors, c("peptide_id", "transition_group_id"))
  }
  end_time <- Sys.time()
  message("The execution time for getting precursors:")
  print(end_time - start_time)
  nativeIDs <- unlist(precursors$transition_ids)
  nativeIDs
}

#' Calculates area of a peak in XIC group
#'
#' Retention time from reference run is mapped to experiment run using AlignObj.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-13
#' @inheritParams checkParams
#' @param XICs (list) list of extracted ion chromatograms of a precursor.
#' @param left (numeric) left boundary of the peak.
#' @param right (numeric) right boundary of the peak.
#' @return area (numeric)
#' @keywords internal
#' @seealso \code{\link{areaIntegrator}, \link{setAlignmentRank}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]]
#' \dontrun{
#' calculateIntensity(XICs, 5220, 5261, integrationType = "intensity_sum",
#'  baselineType = "base_to_base", fitEMG = FALSE)
#' }
calculateIntensity <- function(XICs, left, right, params){
  time <- lapply(XICs, `[`, i =, j =1)
  intensityList <- lapply(XICs, `[`, i =, j= 2)
  if(params[["smoothPeakArea"]]){
    kL <- params[["kernelLen"]]
    pO <- params[["polyOrd"]]
  } else{
    kL <- 0L
    pO <- 1L
  }
  intensity <- areaIntegrator(time, intensityList, left, right,  params[["integrationType"]], params[["baselineType"]],
                              FALSE, params[["baseSubtraction"]], kL, pO)
  intensity[is.nan(intensity)] <- NA_real_
  if(params[["transitionIntensity"]]) return (intensity)
  sum(intensity, na.rm = FALSE)
}


newRow <- function(df, xics, left, right, rt, analyte, Run, params){
  intensity <- calculateIntensity(xics, left, right, params)
  intensity <- ifelse(params[["transitionIntensity"]], list(intensity), intensity)
  idx <- which(df$run == Run & df$transition_group_id == analyte)
  idx <- idx[is.na(.subset2(df, "peak_group_rank")[idx])]
  # idx <- df[run == Run & transition_group_id == analyte, .I[is.na(peak_group_rank)][1], by = run]$V1
  if(length(idx) == 0) return(invisible(NULL))
  set(df, idx[1L], c(3L, 4L, 5L, 6L, 10L), list(rt, intensity, left, right, 1L))
  invisible(NULL)
}


reIntensity <- function(df, Run, XICs, params){
  idx <- df[run == Run & alignment_rank == 1, which = TRUE]
  for(i in idx){
    analyte_chr <- as.character(.subset2(df, "transition_group_id")[[i]])
    area <- calculateIntensity(XICs[[analyte_chr]], .subset2(df, "leftWidth")[[i]], .subset2(df, "rightWidth")[[i]],
                               params)
    data.table::set(df, i, "intensity", area)
  }
  invisible(NULL)
}

#' Calculates area of peaks in peakTable
#'
#' For the give peak boundary in peakTable, the function extracts raw chromatograms and recalculate intensities.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-05-28
#'
#' @import data.table
#' @inheritParams alignTargetedRuns
#' @param peakTable (data-table) usually an output of alignTargetedRuns. Must have these columns: precursor, run, intensity, leftWidth, rightWidth.
#' @return (data-table)
#' @seealso \code{\link{alignTargetedRuns}, \link{calculateIntensity}}
#' @examples
#' library(data.table)
#' peakTable <- data.table(precursor = c(1967L, 1967L, 2474L, 2474L),
#'                    run = rep(c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
#'                    "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"), 2),
#'                    intensity = c(186.166, 579.832, 47.9525, 3.7413),
#'                    leftWidth = c(5001.76, 5025.66, 6441.51, 6516.6),
#'                    rightWidth = c(5076.86, 5121.25, 6475.65, 6554.2))
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' params$smoothPeakArea <- TRUE
#' recalculateIntensity(peakTable, dataPath, params = params)
#' peakTable <- data.table(precursor = c(1967L, 1967L, 2474L, 2474L),
#'                    run = rep(c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
#'                    "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"), 2),
#'                    intensity = list(NA, NA, NA, NA),
#'                    leftWidth = c(5001.76, 5025.66, 6441.51, 6516.6),
#'                    rightWidth = c(5076.86, 5121.25, 6475.65, 6554.2))
#' params$transitionIntensity <- TRUE
#' recalculateIntensity(peakTable, dataPath, params = params)
#' @export
recalculateIntensity <- function(peakTable, dataPath = ".", oswMerged = TRUE, params = paramsDIAlignR()){
  if(!all(colnames(peakTable) == c("precursor", "run", "intensity", "leftWidth", "rightWidth"))){
    stop("Columns must be precursor, run, intensity, leftWidth, rightWidth in the same order.")
  }
  peakTable <- data.table::copy(peakTable)
  data.table::setkeyv(peakTable, "run")
  runs <- unique(peakTable$run)
  analytes <- unique(peakTable$precursor)
  fileInfo <- getRunNames(dataPath, oswMerged, params)
  fileInfo <- updateFileInfo(fileInfo, runs)

  ######### Get Precursors from the query and respectve chromatogram indices. ######
  precursors <- getPrecursorByID(analytes, fileInfo)
  precursors[,peptide_id := NULL]
  data.table::setkeyv(precursors, c("transition_group_id"))

  ######### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(fileInfo)
  message("Metadata is collected from mzML files.")

  ############# Get chromatogram Indices of precursors across all runs. ############
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)
  fetchXICs = extractXIC_group2
  analytes <- precursors[, "transition_group_id"][[1]]
  # Iterate through each run
  for (run in rownames(fileInfo)){
    chromIndices <- prec2chromIndex[[run]][,chromatogramIndex]
    con <- createTemp(mzPntrs[[run]], unlist(chromIndices))
    # Fetch XICs
    XICs <- lapply(seq_along(analytes), function(i){
      cI <- chromIndices[[i]]
      if(any(is.na(unlist(cI))) | is.null(unlist(cI))) return(NULL)
      fetchXICs(con, cI)
    })
    names(XICs) <- as.character(analytes)
    DBI::dbDisconnect(con)
    # Calculate intensity for each analyte
    for(i in seq_along(analytes)){
      idx <- which(peakTable[["run"]] == fileInfo[run, "runName"] & peakTable[["precursor"]] == analytes[i])
      xics <- XICs[[i]]
      if(is.null(xics) || any(vapply(xics, missingInXIC, FALSE, USE.NAMES = FALSE))){
        message("Chromatogram indices for precursor ", analytes[i], " are missing in ", fileInfo[run, "runName"])
        message("Skipping precursor ", analytes[i], " in ", fileInfo[run, "runName"])
        next
      }
      if(length(idx) == 0L) next
      left <- .subset2(peakTable, "leftWidth")[[idx]]
      right <- .subset2(peakTable, "rightWidth")[[idx]]
      intensity <- calculateIntensity(xics, left, right, params)
      intensity <- ifelse(params[["transitionIntensity"]], list(intensity), intensity)
      data.table::set(peakTable, idx, c(3L), intensity)
    }
  }

  for(mz in mzPntrs){
    if(is(mz)[1] == "SQLiteConnection") DBI::dbDisconnect(mz)
    if(is(mz)[1] == "mzRpwiz") rm(mz)
  }
  data.table::setkeyv(peakTable, c("precursor"))
  peakTable
}

reIntensity2 <- function(df, idx, XICs, pk, params){
  area <- calculateIntensity(XICs, pk[1], pk[2], params)
  data.table::set(df, i = idx, c(4L,5L,6L), list(area,pk[1],pk[2]))
}

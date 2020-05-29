#' Calculates area of a peak in XIC group
#'
#' Retention time from reference run is mapped to experiment run using AlignObj.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-13
#'
#' @param XICs (list) list of extracted ion chromatograms of a precursor.
#' @param left (numeric) left boundary of the peak.
#' @param right (numeric) right boundary of the peak.
#' @param integrationType (string) method to ompute the area of a peak contained in XICs. Must be
#'  from "intensity_sum", "trapezoid", "simpson".
#' @param baselineType (string) method to estimate the background of a peak contained in XICs. Must be
#'  from "base_to_base", "vertical_division_min", "vertical_division_max".
#' @param fitEMG (logical) enable/disable exponentially modified gaussian peak model fitting.
#' @param baseSubtraction (logical) TRUE: remove background from peak signal using estimated noise levels.
#' @return (numeric)
#' @keywords internal
#' @seealso \code{\link{getMultipeptide}, \link{setAlignmentRank}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' \dontrun{
#' calculateIntensity(XICs, 5220, 5261, integrationType = "intensity_sum",
#'  baselineType = "base_to_base", fitEMG = FALSE)
#' }
calculateIntensity <- function(XICs, left, right, integrationType, baselineType,
                               fitEMG, baseSubtraction = TRUE){
  time <- lapply(XICs, `[[`, 1)
  intensityList <- lapply(XICs, `[[`, 2)
  intensity <- areaIntegrator(time, intensityList, left, right, integrationType, baselineType,
                              fitEMG, baseSubtraction)
  intensity
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
#' @importFrom magrittr %>%
#' @param peakTable (data-frame) usually an output of alignTargetedRuns. Must have these columns: run, precursor, leftWidth, rightWidth.
#' @param dataPath (string) path to mzml and osw directory.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param XICfilter (string) must be either sgolay, boxcar, gaussian, loess or none.
#' @param polyOrd (integer) order of the polynomial to be fit in the kernel.
#' @param kernelLen (integer) number of data-points to consider in the kernel.
#' @param baseSubtraction (logical) TRUE: remove background from peak signal using estimated noise levels.
#' @param baselineType (string) method to estimate the background of a peak contained in XICs. Must be
#'  from "base_to_base", "vertical_division_min", "vertical_division_max".
#' @param integrationType (string) method to ompute the area of a peak contained in XICs. Must be
#'  from "intensity_sum", "trapezoid", "simpson".
#' @param fitEMG (logical) enable/disable exponentially modified gaussian peak model fitting.
#' @param smoothPeakArea (logical) FALSE: raw chromatograms will be used for quantification. TRUE: smoothed chromatograms will be used for quantification.
#' @return (data-frame)
#' @seealso \code{\link{alignTargetedRuns}, \link{calculateIntensity}}
#' @examples
#' peakTable <- data.frame(precursor = c(1967L, 1967L, 2474L, 2474L),
#'                    run = rep(c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
#'                    "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"), 2),
#'                    intensity = c(186.166, 579.832, 47.9525, 3.7413),
#'                    leftWidth = c(5001.76, 5025.66, 6441.51, 6516.6),
#'                    rightWidth = c(5076.86, 5121.25, 6475.65, 6554.2), stringsAsFactors = FALSE)
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' newTable <- recalculateIntensity(peakTable, dataPath)
#' @export
recalculateIntensity <- function(peakTable, dataPath = ".", oswMerged = TRUE,
                                 XICfilter = "sgolay", polyOrd = 4, kernelLen = 9, baseSubtraction = TRUE,
                                 baselineType = "base_to_base", integrationType = "intensity_sum",
                                 fitEMG = FALSE, smoothPeakArea = FALSE){
  runs <- unique(peakTable$run)
  analytes <- unique(peakTable$precursor)
  fileInfo <- getRunNames(dataPath, oswMerged)
  fileInfo <- updateFileInfo(fileInfo, runs)

  ######### Get Precursors from the query and respectve chromatogram indices. ######
  precursors <- getPrecursorByID(analytes, fileInfo)

  ######### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(fileInfo)
  message("Metadata is collected from mzML files.")

  ############# Get chromatogram Indices of precursors across all runs. ############
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)

  newArea <- list()
  for (run in rownames(fileInfo)){
    newArea[[run]] <- rep(NA_real_, length(analytes))
    runname <- fileInfo[run, "runName"]
    for (i in seq_along(analytes)){
      analyte <- analytes[i]
      df <- dplyr::filter(peakTable, .data$precursor == analyte, .data$run == runname) %>%
        dplyr::select(.data$leftWidth, .data$rightWidth)
      chromIndices <- prec2chromIndex[[run]][["chromatogramIndex"]][[i]]

      # Get XIC_group from reference run. if missing, go to next analyte.
      if(any(is.na(chromIndices))){
        warning("Chromatogram indices for ", analyte, " are missing in ", fileInfo[run, "runName"])
        message("Skipping ", analyte, " in ", fileInfo[run, "runName"], ".")
        next
      } else {
        XICs <- extractXIC_group(mz = mzPntrs[[run]], chromIndices = chromIndices)
      }
      if(smoothPeakArea){
        XICs <- smoothXICs(XICs, type = XICfilter, kernelLen = kernelLen, polyOrd = polyOrd)
      }
      area <- calculateIntensity(XICs, df[1, "leftWidth"], df[1, "rightWidth"],  integrationType = integrationType,
                                 baselineType = baselineType, fitEMG = fitEMG, baseSubtraction = baseSubtraction)
      newArea[[run]][i] <- area
    }
  }

  newArea <- as.data.frame(do.call(cbind, newArea))
  newArea$precursor <- analytes
  newArea <- tidyr::pivot_longer(newArea, -.data$precursor, names_to = "run",
                                 values_to = "intensity") %>% as.data.frame()
  newArea$run <- fileInfo[newArea$run, "runName"]
  newArea
}
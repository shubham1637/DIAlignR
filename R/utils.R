#' Fetch the reference run for each peptide
#'
#' Provides the reference run based on lowest p-value.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-08
#' @importFrom data.table set data.table
#' @inheritParams alignTargetedRuns
#' @param peptideScores (list of data-frames) each dataframe has scores of a peptide across all runs.
#' @return (dataframe) has two columns:
#' \item{peptide_id}{(integer) a unique id for each peptide.}
#' \item{run}{(string) run identifier.}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath)
#' precursorsInfo <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_Proteomics",
#'                                 context = "experiment-wide", maxPeptideFdr = 0.05)
#' peptideIDs <- unique(precursorsInfo$peptide_id)
#' peptidesInfo <- getPeptideScores(fileInfo, peptideIDs)
#' peptidesInfo <- lapply(peptideIDs, function(pep) peptidesInfo [.(pep),])
#' names(peptidesInfo) <- as.character(peptideIDs)
#' getRefRun(peptidesInfo)
#' @seealso \code{\link{getPeptideScores}}
#' @export
getRefRun <- function(peptideScores, applyFun=lapply){
  DFs <- data.table("peptide_id" = as.integer(names(peptideScores)),
             "run" = NA_character_ , key = "peptide_id")
  invisible(
    applyFun(seq_along(peptideScores), function(i){
      DT <- peptideScores[[i]]
      idx <- DT[, which.min(pvalue)]
      if(length(idx)!=0) set(DFs, i, j = 2L, DT[[2]][[idx]])
      invisible(NULL)
    })
  )
  DFs
}

#' Get multipeptides
#'
#' Each element of the multipeptide is a collection of features associated with a peptide.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-08
#' @importFrom data.table rbindlist set setkeyv setcolorder
#' @importFrom bit64 NA_integer64_
#' @inheritParams alignTargetedRuns
#' @param precursors (data-frames) Contains precursors and associated transition IDs.
#' @param features (list of data-frames) Contains features and their properties identified in each run.
#' @param runType (char) This must be one of the strings "DIA_Proteomics", "DIA_IPF", "DIA_Metabolomics".
#' @param masters (characters) names of extra runs.
#' @return (list) of dataframes having following columns:
#' \item{transition_group_id}{(integer) a unique id for each precursor.}
#' \item{run}{(string) run identifier.}
#' \item{RT}{(numeric) retention time as in FEATURE.EXP_RT of osw files.}
#' \item{Intensity}{(numeric) peak intensity as in FEATURE_MS2.AREA_INTENSITY of osw files.}
#' \item{leftWidth}{(numeric) as in FEATURE.LEFT_WIDTH of osw files.}
#' \item{rightWidth}{(numeric) as in FEATURE.RIGHT_WIDTH of osw files.}
#' \item{peak_group_rank}{(integer) rank of each feature associated with transition_group_id.}
#' \item{m_score}{(numeric) q-value of each feature associated with transition_group_id.}
#' \item{alignment_rank}{(integer) rank of each feature post-alignment.}
#'
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
#' precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide")
#' features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
#' multipeptide <- getMultipeptide(precursors, features)
#' multipeptide[["9861"]]
#' @seealso \code{\link{getPrecursors}, \link{getFeatures}}
#' @export
getMultipeptide <- function(precursors, features, runType="DIA_Proteomics", applyFun=lapply, masters = NULL){
  peptideIDs <- unique(precursors$peptide_id)
  runs <- names(features)
  num_run = length(features)
  if(is.null(masters)){
    multipeptide <- applyFun(seq_along(peptideIDs), function(i){
      # Get transition_group_id for a peptide
      analytes <- precursors[list(peptideIDs[i]), 1L][[1]]
      num_analytes <- length(analytes)
      newdf <- rbindlist(lapply(runs, function(run){
        df <- rbindlist(list(features[[run]][list(analytes), ],
                  dummyTbl(analytes, runType)), use.names=TRUE) # dummy for new features
        set(df, j = c("run", "alignment_rank"), value = list(run, NA_integer_))
        df
      }), use.names=TRUE)
      if ( runType=="DIA_IPF" ){
        setcolorder(newdf, c('transition_group_id', 'feature_id', 'RT', 'intensity', 'leftWidth', 'rightWidth', 'peak_group_rank', 'm_score', 'run', 'alignment_rank', 'ms2_m_score'))
      }
      setkeyv(newdf, "run")
      newdf
    })
  } else{
    multipeptide <- applyFun(seq_along(peptideIDs), function(i){
      # Get transition_group_id for a peptide
      analytes <- precursors[list(peptideIDs[i]), 1L][[1]]
      num_analytes <- length(analytes)
      df1 <- lapply(runs, function(run){
        df <- features[[run]][list(analytes), ]
        df[,`:=`("run" = run, "alignment_rank" = NA_integer_)]
        df
      })
      df2 <- dummyMerge(analytes, masters, runType)
      newdf <- rbindlist(list(rbindlist(df1, use.names=TRUE), df2), use.names=TRUE)
      if ( runType=="DIA_IPF" ){
        setcolorder(newdf, c('transition_group_id', 'feature_id', 'RT', 'intensity', 'leftWidth', 'rightWidth', 'peak_group_rank', 'm_score', 'run', 'alignment_rank', 'ms2_m_score'))
      }
      setkeyv(newdf, "run")
      newdf
    })
  }

  # Convert peptides as character. Add names to the multipeptide list.
  names(multipeptide) <- as.character(peptideIDs)
  multipeptide
}

dummyTbl <- function(analytes, runType="DIA_Proteomics"){
  if ( runType=="DIA_IPF" ){
    data.table("transition_group_id" = analytes, "feature_id" = bit64::NA_integer64_,
               "RT" = NA_real_, "intensity" = NA_real_, "leftWidth" = NA_real_,
               "rightWidth" = NA_real_, "peak_group_rank" = NA_integer_, "ms2_m_score" = NA_real_, "m_score" = NA_real_)
  } else {
    data.table("transition_group_id" = analytes, "feature_id" = bit64::NA_integer64_,
               "RT" = NA_real_, "intensity" = NA_real_, "leftWidth" = NA_real_,
               "rightWidth" = NA_real_, "peak_group_rank" = NA_integer_, "m_score" = NA_real_)
  }

}

dummyMerge <- function(analytes, masters, runType="DIA_Proteomics"){
  if ( runType=="DIA_IPF" ){
    data.table("transition_group_id" = rep(analytes, times = length(masters), each = 5L), "feature_id" = bit64::NA_integer64_,
               "RT" = NA_real_, "intensity" = NA_real_, "leftWidth" = NA_real_,
               "rightWidth" = NA_real_, "peak_group_rank" = NA_integer_, "ms2_m_score" = NA_real_, "m_score" = NA_real_,
               "run" = rep(masters, each = 5L*length(analytes)), "alignment_rank" = NA_integer_)
  } else {
    data.table("transition_group_id" = rep(analytes, times = length(masters), each = 5L), "feature_id" = bit64::NA_integer64_,
               "RT" = NA_real_, "intensity" = NA_real_, "leftWidth" = NA_real_,
               "rightWidth" = NA_real_, "peak_group_rank" = NA_integer_, "m_score" = NA_real_,
               "run" = rep(masters, each = 5L*length(analytes)), "alignment_rank" = NA_integer_)
  }
}

#' Writes the output table post-alignment
#'
#' Selects all features from multipeptide with alignment rank = 1, and write them onto the disk.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-14
#' @importFrom data.table rbindlist setnames setorder setcolorder set
#' @param filename (string) Name of the output file.
#' @param fileInfo (data-frame) Output of getRunNames function.
#' @param multipeptide (list of data-frames) Each element of the list is collection of features associated with a precursor.
#' @param precursors (data-frame) for each transition_group_id, contains peptide sequence and charge.
#' @return An output table with following columns: precursor, run, intensity, RT, leftWidth, rightWidth,
#'  peak_group_rank, m_score, alignment_rank, peptide_id, sequence, charge, group_label.
#' @seealso \code{\link{getRunNames}, \link{getMultipeptide}, \link{getPrecursors}}
#' @keywords internal
#'
#' @examples
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' \dontrun{
#' writeTables(fileInfo, multipeptide, precursors)
#' }
writeTables <- function(fileInfo, multipeptide, precursors){
  peptides <- precursors[, logical(1), keyby = peptide_id]$peptide_id
  runs <- rownames(fileInfo)
  idx <- grep("^master[0-9]+$", runs, invert = TRUE)
  runs <- runs[idx]
  runName <- fileInfo[idx, "runName"]

  #### Get a dataframe of all analytes with alignment rank = 1 ###########
  finalTbl <- lapply(seq_along(peptides), function(i){
    indices <- rep(NA_integer_, length(runs))
    df <- multipeptide[[i]]
    idx <- unlist(lapply(runs, function(run){
      j <- which(df[["run"]] == run)
      idx <- j[df$alignment_rank[j]  == 1L]
      #TODO: If alignment-rank is 1 but out of boundary then select other peak?
      if(all(is.na(idx))) idx <- j[df$peak_group_rank[j] == 1L]
      idx
    }))
    idx <- idx[!is.na(idx)]
    df[idx,]
  })
  finalTbl <- rbindlist(finalTbl)
  finalTbl[,run := runName[match(finalTbl$run, runs)]]

  ##### Merging precursor information and return the dataframe. ###################
  finalTbl <- precursors[, !("transition_ids")][finalTbl, on = "transition_group_id"]
  setnames(finalTbl, "transition_group_id", "precursor")
  set(finalTbl, j = "feature_id", value = as.character(finalTbl[["feature_id"]]))

  setorder(finalTbl, peptide_id, precursor, run)
  setcolorder(finalTbl, c("peptide_id", "precursor", "run", "RT", "intensity", "leftWidth", "rightWidth",
                    "peak_group_rank", "m_score", "alignment_rank", "feature_id"))
  finalTbl
}

# Alignment of precursor 4618 , sequence = QFNNTDIVLLEDFQK/3 across runs
# ref = hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt
# eXp = hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"
# Example: test_getAlignObj
testAlignObj <- function(){
  AlignObj <- new("AffineAlignObjLight",
                    indexA_aligned = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,0,20,0,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,0,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,0,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,0,162,0,163,164,165,166,0,167,0,168,169,170,171,172,173,174,175,176),
                    indexB_aligned = c(0,0,0,1,0,2,0,3,0,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,0,92,0,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176),
                    score = c(0,0,0,2.675751,2.385165,4.745081,4.454496,6.67706,6.386474,9.641135,15.48268,26.87968,46.67034,77.78939,121.8607,170.8698,214.2358,244.4554,262.8537,262.5631,270.8538,270.5632,288.2477,319.3287,364.6418,413.4873,453.1931,479.6588,496.5207,506.6607,512.7442,515.267,516.8242,518.2747,519.8424,521.2872,522.6472,524.2912,525.6285,526.5892,526.9713,527.1531,527.4022,527.1116,530.9457,547.7525,588.2834,658.8819,748.3079,833.8337,898.3289,935.5809,948.8015,952.0709,952.8035,953.4267,954.0863,954.8143,955.4842,956.0834,956.802,957.535,958.2853,959.0355,959.7972,960.7983,961.8922,963.0142,964.2597,965.5837,966.878,968.0037,968.4412,968.1507,968.1958,968.9242,985.3144,1085.833,1364.976,1846.928,2409.31,2869.416,3132.509,3231.061,3257.015,3264.422,3269.377,3275.003,3282.515,3290.524,3297.864,3304.43,3310.324,3314.403,3316.806,3317.992,3318.933,3318.642,3319.328,3319.038,3320.17,3321.781,3323.71,3325.64,3327.855,3330.382,3332.989,3335.319,3337.555,3339.96,3342.381,3344.48,3346.456,3348.605,3350.446,3352.092,3353.829,3355.911,3358.256,3360.576,3363.292,3367.099,3372.687,3380.124,3389.957,3401.498,3414.81,3428.762,3441.046,3451.052,3459.235,3466.392,3473.212,3480.14,3490.173,3506.584,3530.062,3561.003,3595.718,3624.828,3642.574,3650.352,3653.893,3656.295,3658.798,3661.361,3663.704,3665.936,3667.714,3669.478,3670.721,3671.991,3673.278,3674.689,3676.068,3677.317,3678.688,3680.062,3681.513,3683.097,3684.786,3686.565,3688.24,3689.741,3690.859,3690.568,3691.496,3691.205,3692.31,3693.465,3694.458,3695.352,3695.061,3695.892,3695.602,3696.512,3697.468,3698.18,3698.799,3699.363,3699.94,3700.634,3701.585,3702.988))

  AlignObj
}

#' Checks all the alignment parameters
#'
#' Parameters are defined in the \code{\link{paramsDIAlignR}} function. If parameters are found
#' incongruous, the program will terminate.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-01
#' @param params (list) parameters are entered as list. Output of the \code{\link{paramsDIAlignR}} function.
#' @return Invisible NULL
#' @examples
#' params <- paramsDIAlignR()
#' \dontrun{
#' checkParams(params)
#' }
#' @seealso \code{\link{paramsDIAlignR}}
#' @keywords internal
checkParams <- function(params){
  if(!(params[["chromFile"]] %in% c("mzML", "sqMass"))) stop("chromFile must either be mzML or sqMass.")

  if(params[["context"]] != "experiment-wide" & params[["context"]] != "global"){
    stop("context must either be experiment-wide or global for the alignment.")
  }

  if(params[["maxFdrQuery"]]>1 | params[["maxFdrQuery"]] < 0){
    stop("maxFdrQuery must be between 0 and 1. Recommended value is 0.05")
  }

  if(params[["maxIPFFdrQuery"]]>1 | params[["maxIPFFdrQuery"]] < 0){
    stop("maxIPFFdrQuery must be between 0 and 1. Recommended value is 0.05")
  }

  if(params[["maxPeptideFdr"]] < 0 | params[["maxPeptideFdr"]] > 1){
    stop("maxPeptideFdr must be between 0 and 1. Recommended value is 0.01")
  }

  if(!any(params[["XICfilter"]] %in% c("sgolay",  "boxcar", "gaussian", "loess", "none"))){
    stop("XICfilter must either be sgolay, boxcar, gaussian, loess or none.")
  }

  if(params[["XICfilter"]] == "sgolay"){
    if( (params[["kernelLen"]] %% 2) != 1){
      stop("For XICfilter = sgolay, kernelLen can only be odd number.")
    }
    if(params[["polyOrd"]] > params[["kernelLen"]]) {
      stop("For XICfilter = sgolay, polyOrd must be less than kernelLen.")
    }
  }

  if(params[["XICfilter"]] == "none") params[["kernelLen"]] <- 0L

  if(!any(params[["globalAlignment"]] %in% c("loess",  "linear"))){
    stop("globalAlignment must either be loess or linear.")
  }

  if(params[["globalAlignment"]] == "loess" & (params[["globalAlignmentSpan"]] >1 | params[["globalAlignmentSpan"]] < 0)){
    stop("For globalAlignment = loess, globalAlignmentSpan must be between 0 and 1.")
  }

  if(params[["globalAlignmentFdr"]] < 0 | params[["globalAlignmentFdr"]] > params[["maxFdrQuery"]]){
    stop("globalAlignmentFdr must be between 0 and maxFdrQuery.")
  }

  if(!any(params[["normalization"]] %in% c("mean", "l2", "none"))){
    stop("normalization must be either mean, l2 or none.")
  }

  if(!any(params[["simMeasure"]] %in% c("dotProduct", "cosineAngle", "crossCorrelation", "cosine2Angle",
                            "dotProductMasked", "euclideanDist", "covariance", "correlation"))){
    stop("Choose appropriate simMeasure.")
  }

  if(!any(params[["alignType"]] %in% c("global", "local", "hybrid"))){
    stop("alignType must be either global, local or hybrid. Recommended value is hybrid.")
  }

  if(params[["unalignedFDR"]] < 0 | params[["unalignedFDR"]] > params[["maxFdrQuery"]]){
    stop("unalignedFDR must be between 0 and maxFdrQuery. Recommended value is 0.01")
  }

  if(params[["alignedFDR2"]] < params[["unalignedFDR"]] | params[["alignedFDR2"]] > params[["maxFdrQuery"]]){
    stop("alignedFDR2 must be between unalignedFDR and maxFdrQuery. Recommended value is 0.05")
  }

  if(!any(params[["integrationType"]] %in% c("intensity_sum", "trapezoid", "simpson"))){
    stop("integrationType must be either intensity_sum, trapezoid or simpson.")
  }

  if(!any(params[["baselineType"]] %in% c("none","base_to_base", "vertical_division_min", "vertical_division_max"))){
    stop("baselineType must be either base_to_base, vertical_division_min or vertical_division_max.")
  }

  # DOI: 10.1002/cem.1343
  if(params[["fitEMG"]] == TRUE){
    warning("Exponentially Modified Gaussian is not implemented yet.")
  }

  if(params[["analyteFDR"]] < 0 | params[["analyteFDR"]] > 1){
    # Not used yet
  }

  if(params[["fractionNum"]] < 1){
    stop("Number of fractions must be greater than 1.")
  }

  if(params[["hardConstrain"]]){
    params[["samples4gradient"]] <- 1L
  }

  if(params[["fraction"]] < 1 | params[["fraction"]] > params[["fractionNum"]]){
    stop("fraction must be between 1 and fractionNum.")
  }

  if(!any(params[["level"]] %in% c("Peptide", "Protein"))){
    stop("level for maxPeptideFDR should be either Peptide or Protein.")
  }

  if(params[["useIdentifying"]] != TRUE & params[["useIdentifying"]] != FALSE){
    stop("useIdentifying can only be TRUE or FALSE.")
  }

  params
}

#' Parameters for the alignment functions
#'
#' Retention alignment requires OpenSWATH/pyProphet extracted features and chromatograms. This function provides
#' a suite of parameters used for selecting features and manipulating chromatograms. Chromatogram
#' alignment can be performed via reference based or progressively via rooted or unrooted tree. This
#' function provides sensible parameters for these tasks.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-11
#' @return A list of parameters:
#' \item{runType}{(string) must be one of the strings "DIA_Proteomics", "DIA_IPF", "DIA_Metabolomics".}
#' \item{chromFile}{(string) must either be "mzML" or "sqMass".}
#' \item{maxFdrQuery}{(numeric) a numeric value between 0 and 1. It is used to filter peptides from osw file which have SCORE_MS2.QVALUE less than itself.}
#' \item{maxIPFFdrQuery}{(numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_IPF.QVALUE less than itself. (For PTM IPF use)}
#' \item{maxPeptideFdr}{(numeric) a numeric value between 0 and 1. It is used to filter peptides from osw file which have SCORE_PEPTIDE.QVALUE less than itself.}
#' \item{analyteFDR}{(numeric) the upper limit of feature FDR to be it considered for building tree.}
#' \item{treeDist}{(string) the method used to build distance matrix. Must be either "rsquared" or "count.}
#' \item{context}{(string) used in pyprophet peptide. Must be either "run-specific", "experiment-wide", or "global".}
#' \item{unalignedFDR}{(numeric) must be between 0 and maxFdrQuery. Features below unalignedFDR are
#'  considered for quantification even without the RT alignment.}
#' \item{alignedFDR1}{(numeric) must be between unalignedFDR and alignedFDR2. Features below alignedFDR1 and aligned to the reference are
#'  considered for quantification.}
#' \item{alignedFDR2}{(numeric) must be between alignedFDR1 and maxFdrQuery. Features below alignedFDR2 and within certain distance from the aligned time are
#'  considered for quantification after the alignment.}
#' \item{level}{(string) apply maxPeptideFDR on Protein as well if specified as "Protein". Default: "Peptide".}
#' \item{integrationType}{(string) method to ompute the area of a peak contained in XICs. Must be
#'  from "intensity_sum", "trapezoid", "simpson".}
#' \item{baseSubtraction}{{logical} TRUE: remove background from peak signal using estimated noise levels.}
#' \item{baselineType}{(string) method to estimate the background of a peak contained in XICs. Must be
#'  from "none", "base_to_base", "vertical_division_min", "vertical_division_max".}
#' \item{fitEMG}{(logical) enable/disable exponentially modified gaussian peak model fitting.}
#' \item{recalIntensity}{(logical) recalculate intensity for all analytes.}
#' \item{fillMissing}{(logical) calculate intensity for ananlytes for which features are not found.}
#' \item{XICfilter}{(string) must be either sgolay, boxcar, gaussian, loess or none.}
#' \item{polyOrd}{(integer) order of the polynomial to be fit in the kernel.}
#' \item{kernelLen}{(integer) number of data-points to consider in the kernel.}
#' \item{globalAlignment}{(string) must be either "loess" or "linear".}
#' \item{globalAlignmentFdr}{(numeric) a numeric value between 0 and 1. Features should have m-score lower than this value for participation in LOESS fit.}
#' \item{globalAlignmentSpan}{(numeric) spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.}
#' \item{RSEdistFactor}{(numeric) defines how much distance in the unit of rse remains a noBeef zone.}
#' \item{normalization}{(string) must be selected from "mean", "l2".}
#' \item{simMeasure}{(string) must be selected from dotProduct, cosineAngle, crossCorrelation,
#'   cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.}
#' \item{alignType}{(numeric) available alignment methods are "global", "local" and "hybrid".}
#' \item{goFactor}{(numeric) penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.}
#' \item{geFactor}{(numeric) penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.}
#' \item{cosAngleThresh}{(numeric) in simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.}
#' \item{OverlapAlignment}{(logical) an input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.}
#' \item{dotProdThresh}{(numeric) in simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.}
#' \item{gapQuantile}{(numeric) must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.}
#' \item{kerLen}{(integer) In simType = crossCorrelation, length of the kernel used to sum similarity score. Must be an odd number.}
#' \item{hardConstrain}{(logical) if FALSE; indices farther from noBeef distance are filled with distance from linear fit line.}
#' \item{samples4gradient}{(numeric) modulates penalization of masked indices.}
#' \item{fillMethod}{(string) must be either "spline", "sgolay" or "linear".}
#' \item{splineMethod}{(string) must be either "fmm" or "natural".}
#' \item{mergeTime}{(string) must be either "ref", "avg", "refStart" or "refEnd".}
#' \item{keepFlanks}{(logical) TRUE: Flanking chromatogram is not removed.}
#' \item{fraction}{(integer) indicates which fraction to align.}
#' \item{fractionNum}{(integer) Number of fractions to divide the alignment.}
#' \item{lossy}{(logical) if TRUE, time and intensity are lossy-compressed in generated sqMass file.}
#' \item{useIdentifying}{(logical) Set TRUE to use identifying transitions in alignment. (DEFAULT: FALSE)}
#' @seealso \code{\link{checkParams}, \link{alignTargetedRuns}}
#' @examples
#' params <- paramsDIAlignR()
#' @export
paramsDIAlignR <- function(){
  params <- list( runType = "DIA_Proteomics", chromFile = "sqMass",
                  maxFdrQuery = 0.05, maxIPFFdrQuery = 0.05, maxPeptideFdr = 0.01, analyteFDR = 0.01, treeDist = "rsquared",
                  context = "global", unalignedFDR = 0.01, alignedFDR1 = 0.01, alignedFDR2 = 0.05, level = "Peptide",
                  integrationType = "intensity_sum", baselineType = "base_to_base", fitEMG = FALSE,
                  recalIntensity = FALSE, fillMissing = TRUE, baseSubtraction = TRUE,
                  XICfilter = "sgolay", polyOrd = 4L, kernelLen = 11L,
                  globalAlignment = "loess", globalAlignmentFdr = 0.01, globalAlignmentSpan = 0.1,
                  RSEdistFactor = 3.5, normalization = "mean", simMeasure = "dotProductMasked",
                  alignType = "hybrid", goFactor = 0.125, geFactor = 40,
                  cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                  dotProdThresh = 0.96, gapQuantile = 0.5, kerLen = 9,
                  hardConstrain = FALSE, samples4gradient = 1L,
                  fillMethod = "spline", splineMethod = "natural", mergeTime = "avg", smoothPeakArea = FALSE,
                  keepFlanks = TRUE, wRef = 0.5, batchSize = 1000L, transitionIntensity = FALSE,
                  fraction = 1L, fractionNum = 1L, lossy = FALSE, useIdentifying = FALSE)
  params
}

#' Prints alignment summary
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-15
#' @param finalTbl (dataframe) an output of  \code{\link{writeTables}} function.
#' @param params (list) must have following elements: alignedFDR2 and unalignedFDR.
#' @return Invisible NULL
#' @seealso \code{\link{paramsDIAlignR}, \link{writeTables}}
#' @keywords internal
alignmentStats <- function(finalTbl, params){
  if (params[["runType"]]=="DIA_IPF"){
    # Without alignment at unaligned FDR:
    num1 <- sum(finalTbl$original_m_score <= params[["unalignedFDR"]] & finalTbl$alignment_rank == 1L &
                  !is.na(finalTbl$intensity), na.rm = TRUE)
    message("The number of quantified precursors at ", params[["unalignedFDR"]], " FDR: ", num1)

    # Without alignment at aligned FDR (Gain):
    num2 <- sum(finalTbl$original_m_score > params[["unalignedFDR"]] & finalTbl$alignment_rank == 1L &
                  !is.na(finalTbl$intensity), na.rm = TRUE)
    message("The increment in the number of quantified precursors at ", params[["alignedFDR2"]],
            " FDR: ", num2)

    # Corrected peptides by Alignmet
    idx <- finalTbl$peak_group_rank != 1L & finalTbl$original_m_score > params[["unalignedFDR"]] &
      finalTbl$alignment_rank == 1L & !is.na(finalTbl$intensity)
    num3 <- sum(idx, na.rm = TRUE)
    message("Out of ", num2, " DIAlignR corrects the peaks for ", num3, " precursors.")
    message("Hence, it has corrected quantification of ", round(num3*100/(num2 + num1), 3), "% precursors.")

    # Gain by calculating area of missing features:
    num4 <- sum(!is.na(finalTbl$intensity) & is.na(finalTbl$original_m_score) & finalTbl$alignment_rank == 1L, na.rm = TRUE)
    message("DIAlignR has calculated quantification for ", num4, " precursors, for which peaks were not identified.")
    message("Thus, it provides a gain of ", round(num4*100/(num2 + num1 + num4), 3), "%.")

    if(params[["fillMissing"]]) message(sum(is.na(finalTbl$intensity)),
                                        " precursors had part of the aligned peak out of the chromatograms or missing chromatograms, hence could not be quantified.")

  } else {
    # Without alignment at unaligned FDR:
    num1 <- sum(finalTbl$m_score <= params[["unalignedFDR"]] & finalTbl$alignment_rank == 1L &
                  !is.na(finalTbl$intensity), na.rm = TRUE)
    message("The number of quantified precursors at ", params[["unalignedFDR"]], " FDR: ", num1)

    # Without alignment at aligned FDR (Gain):
    num2 <- sum(finalTbl$m_score > params[["unalignedFDR"]] & finalTbl$alignment_rank == 1L &
                  !is.na(finalTbl$intensity), na.rm = TRUE)
    message("The increment in the number of quantified precursors at ", params[["alignedFDR2"]],
            " FDR: ", num2)

    # Corrected peptides by Alignmet
    idx <- finalTbl$peak_group_rank != 1L & finalTbl$m_score > params[["unalignedFDR"]] &
      finalTbl$alignment_rank == 1L & !is.na(finalTbl$intensity)
    num3 <- sum(idx, na.rm = TRUE)
    message("Out of ", num2, " DIAlignR corrects the peaks for ", num3, " precursors.")
    message("Hence, it has corrected quantification of ", round(num3*100/(num2 + num1), 3), "% precursors.")

    # Gain by calculating area of missing features:
    num4 <- sum(!is.na(finalTbl$intensity) & is.na(finalTbl$m_score) & finalTbl$alignment_rank == 1L, na.rm = TRUE)
    message("DIAlignR has calculated quantification for ", num4, " precursors, for which peaks were not identified.")
    message("Thus, it provides a gain of ", round(num4*100/(num2 + num1 + num4), 3), "%.")

    if(params[["fillMissing"]]) message(sum(is.na(finalTbl$intensity)),
                            " precursors had part of the aligned peak out of the chromatograms or missing chromatograms, hence could not be quantified.")
  }
  invisible(NULL)
}

#' Prints messages if a certain number of analytes are aligned
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-26
#' @return Invisible NULL
#' @keywords internal
updateOnalignTargetedRuns <-  function(i){
  if(i < 5){
    message(i, " peptides have been aligned.")
  } else if(i < 1000){
    if(i %% 100 == 0) message(i, " peptides have been aligned.")
  } else {
    if(i %% 1000 == 0) message(i, " peptides have been aligned.")
  }
  invisible(NULL)
}

envName <- "TricEnvr"
# helper function to skip tests if we don't have the 'foo' module
skip_if_no_pyopenms <- function() {
  ropenms <- try(get_ropenms(condaEnv = envName, useConda=TRUE))
  no_ropenms <- is(ropenms, "try-error")
  if(no_ropenms)
    testthat::skip("ropenms not available for testing. A conda environment with name TricEnvr is MUST for testing.")
}


#' Overlap of two time ranges
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-17
#' @param x (numeric) must have only two values and x[2] > x[1].
#' @param y (numeric) must have only two values and y[2] > y[1].
#' @return (logical) TRUE: both time ranges overlap. FALSE: both time ranges do not overlap.
#' @keywords internal
#' @seealso getChildFeature
#' @examples
#' \dontrun{
#' checkOverlap(c(9.1, 13.1), c(2.1, 3.1))
#' checkOverlap(c(1.1, 3.1), c(3.2, 7.1))
#' }
checkOverlap <- function(x, y){
  if(any(is.na(c(x,y)))) return(NA)
  # y has left boundary between x
  leftOverlap <- (y[1] - x[1]) >= 0 & (y[1] - x[2]) <= 0
  # y has right boundary between x
  rightOverlap <- (y[2] - x[1]) >= 0 & (y[2] - x[2]) <= 0
  # y is over-arching x
  overArch <- (y[2] - x[2]) >= 0 & (y[1] - x[1]) <= 0
  olap <- (leftOverlap | rightOverlap | overArch)
  olap
}

overlapLen <- function(df, pk, idx){
  left <- pmax(pk[1], .subset2(df, "leftWidth")[idx])
  right <- pmin(pk[2], .subset2(df, "rightWidth")[idx])
  right - left
}


#' Prints messages if a certain number of analytes are aligned
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-11-19
#' @return Invisible NULL
#' @keywords internal
getPrecursorSubset <- function(precursors, params){
  peptideIDs <- unique(precursors$peptide_id)
  len = length(peptideIDs)
  a = params[["fraction"]]
  b = params[["fractionNum"]]
  pepStart <- floor((a-1)*len/b)+1
  pepEnd <- min(floor(a*len/b), len)
  r <- rle(precursors$peptide_id)
  if(a == 1) {
    pepStart = 1
  } else{
    pepStart <- sum(r$lengths[1:(which(r$value == peptideIDs[pepStart])-1)], 1)
  }
  pepEnd <- sum(r$lengths[1:which(r$value== peptideIDs[pepEnd])])
  c(pepStart, pepEnd)
}

#' Re-Assign FDR Aligned Peaks
#'
#' This method will re-assign peaks with either the MS2 FDR if there was not enough confidence (identifying transitions) in the IPF FDR, or it will assign the user defined alignedFDR2. This assumes the reference run has a really good peak that meets the IPF FDR threshold.
#'
#' @author Justin Sing, \email{justinc.sing@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-0386-0092
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-17
#' @param dt (data.table) Aligned results tables from writeTables
#' @param refRuns (dataframe) A peptide reference run table from getRefRun
#' @param fileInfo (dataframe) A table with run filename information from getRunNames
#' @param params (list) A list of parameters used by DIAlignR generated from paramsDIAlignR
#' @return (data.table) Aligned results table with reassigned m-scores
#' @keywords internal
#' @seealso writeTables, getRefRun, getRunNames, paramsDIAlignR
#' @import data.table
#' @importFrom data.table merge.data.table setnames fifelse ':='
#' @examples
#' \dontrun{
#' finalTbl <- ipfReassignFDR(finalTbl, refRuns, fileInfo, params)
#' }
ipfReassignFDR <- function(dt, refRuns, fileInfo, params){
  ## Create Reference Run table with runName to map back to dt
  runNameMap <- fileInfo$runName
  names(runNameMap) <- rownames(fileInfo)
  refRuns[, ref_run := run, by=1:nrow(refRuns)]
  ## Remap from run number id to experiment name id
  refRuns[ .(ref_run = names(runNameMap), to = runNameMap), on = "ref_run", ref_run := i.to]
  ## Remove run numnber id
  refRuns[, run:=NULL]
  ## Merge peptide Reference Run table with aligned results table
  dt <- merge(dt, refRuns, by="peptide_id")
  ## Assign new FDRs
  ## 1. For aligned peaks that had a poor IPF FDR, assign MS2 FDR
  ## 2. For aligned peaks that have a poor IPF FDR and poor MS2 FDR, assign user defined alignedFDR2
  ## TODO: Is this it reasonable to reassign FDRs like this?
  dt[, m_score_new := fifelse((ms2_m_score < m_score & run!=ref_run), ms2_m_score, m_score, na=params[["alignedFDR2"]]), by = 1:nrow(dt)]
  dt[, m_score_new := fifelse((m_score_new < params[["alignedFDR2"]]), m_score_new, params[["alignedFDR2"]], na=params[["alignedFDR2"]]), by = 1:nrow(dt)]
  setnames(dt, c("m_score_new", "m_score"), c("m_score", "original_m_score"))
  dt
}

missingInXIC <- function(XICs){
  any(sapply(seq_along(XICs), function(i) any(is.na(XICs[[i]]))))
}

distMatrix <- function(features, params, applyFun = lapply){
  if(params[["treeDist"]] == "rsquared"){
    distMat <- distMat.RT(features, params, applyFun)
  } else{
    distMat <- distMat.count(features, params, applyFun)
  }
  stats::as.dist(distMat, diag = FALSE, upper = FALSE)
}

distMat.RT <- function(features, params, applyFun = lapply){
  runs <- names(features)
  distMat <- sapply(runs, function(ref){
    c(applyFun(runs, function(eXp){
      RUNS_RT <- tryCatch(expr = getRTdf(features, ref, eXp, params[["analyteFDR"]]),
                          error = function(c) matrix(c(0,0), nrow=1L)
      )
      n <- nrow(RUNS_RT)
      d <- 1.0
      if(n > 2L){
        r <- stats::cor(RUNS_RT[,RT.ref], RUNS_RT[,RT.eXp], method = "pearson")
        d <- min(1, (1-r*r)*(n-1)/(n-2))
      }
      d
    }))
  })
  as.matrix(distMat)
}

distMat.count <- function(features, params, applyFun = lapply){
  ##### Get distances among runs based on the number of high-quality features. #####
  tmp <- applyFun(features, function(df)
    df[df[["m_score"]] <= params[["analyteFDR"]] & df[["peak_group_rank"]] == 1L, "transition_group_id"][[1]])
  #tmp <- tmp[order(names(tmp), decreasing = FALSE)]
  allIDs <- unique(unlist(tmp, recursive = FALSE, use.names = TRUE))
  allIDs <- sort(allIDs)
  simMat <- crossprod(table(utils::stack(tmp))) # Number of common peaks in each pair
  distMat <- length(allIDs) - simMat
  distMat <- distMat/length(allIDs)
  #distMat <- stats::dist(distMat, method = "manhattan")
  distMat
}

#' Inverted versions of in, is.null and is.na
#'
#' @noRd
#'
#' @examples
#' 1 %not_in% 1:10
#' not_null(NULL)
`%not_in%` <- Negate(`%in%`)

not_null <- Negate(is.null)

not_na <- Negate(is.na)

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

#' Alignment Feature Mapping Table
#'
#' Construct an alignment feature mapping table for star align method, to map aligned features in experiment runs to the reference run
#' @author Justin Sing, \email{justinc.sing@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-0386-0092
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2022-11-07
#' @import data.table
#' @inheritParams alignTargetedRuns
#' @inheritParams getMultipeptide
#' @param features (list of data-frames) Contains features and their properties identified in each run.
#' @return (list) of dataframes per precursor peptide id having following columns:
#' \item{reference_feature_id}{(integer) id of reference feature.}
#' \item{experiment_feature_id}{(integer) id of experiment feature aligned to reference feature.}
#'
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
#' precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide")
#' features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
#' multiFeatureAlignmentMap <- getRefExpFeatureMap(precursors, features)
#' multiFeatureAlignmentMap[["9861"]]
#' @seealso \code{\link{getPrecursors}, \link{getFeatures}}
#' @export
getRefExpFeatureMap <- function(precursors, features, applyFun=lapply){
  peptideIDs <- unique(precursors$peptide_id)
  runs <- names(features)
  num_run = length(features)
  # Generate mapping tables per peptide
  multiFeatureAlignmentMap <- lapply(seq_along(peptideIDs), function(i){
    analytes <- precursors[list(peptideIDs[i]), 1L][[1]]
    # Initiate empty data.table with max number of rows to fill
    # TODO: Assume 5 peak groups features per peptide precursor for now. This is the default number of features reported by OpenSwathWorkflow
    n_top_features <- 5
    n_rows <- (length(analytes) * n_top_features) * num_run
    feature_alignment_mapping <- data.table(reference_feature_id=rep(0, n_rows), experiment_feature_id=rep(0, n_rows))
    # Enforce class for columns to be interger64
    data.table::setattr(feature_alignment_mapping[['reference_feature_id']], "class","integer64")
    data.table::setattr(feature_alignment_mapping[['experiment_feature_id']], "class","integer64")
    data.table::setkeyv(feature_alignment_mapping, "reference_feature_id")
    feature_alignment_mapping
  })
  names(multiFeatureAlignmentMap) <- as.character(peptideIDs)
  multiFeatureAlignmentMap
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
  idx <- grep("^run[0-9]+$", runs)
  # idx <- grep("^master[A-Za-z0-9]+$", runs, invert = TRUE)
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

#' Write out alignment map to disk
#'
#' Save alignment mapping to disk, either append table to OSW file, or save TSV file(s)
#' @author Justin Sing, \email{justinc.sing@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-0386-0092
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2022-11-07
#' @import data.table RSQLite
#' @inheritParams alignTargetedRuns
#' @inheritParams getMZMLpointers
#' @param multiFeatureAlignmentMap (list) contains multiple data-frames that are collection of experiment feature ids
#' mapped to corresponding reference feature id per analyte. This is an output of \code{\link{getRefExpFeatureMap}}.
#' @export
writeOutFeatureAlignmentMap <- function(multiFeatureAlignmentMap, oswMerged, fileInfo)
{
  RefExpFeatureMap <- data.table::rbindlist(multiFeatureAlignmentMap)
  RefExpFeatureMap <- RefExpFeatureMap[reference_feature_id!=0 | experiment_feature_id!=0]

  data.table::setorder(RefExpFeatureMap, "reference_feature_id")
  RefExpFeatureMap[, ALIGNMENT_GROUP_ID := .GRP, by = "reference_feature_id"]

  RefExpFeatureMap = data.table::melt(RefExpFeatureMap, if.col = "ALIGNMENT_GROUP_ID", measure.vars = c("reference_feature_id", "experiment_feature_id"), variable.name="REFERENCE", value.name = "FEATURE_ID" )
  RefExpFeatureMap[, REFERENCE:=as.character(REFERENCE)]
  RefExpFeatureMap[.(REFERENCE = c("reference_feature_id", "experiment_feature_id"), to = c("1", "0")), on = "REFERENCE", REFERENCE := i.to]
  RefExpFeatureMap[, ALIGNMENT_GROUP_ID:=as.integer(ALIGNMENT_GROUP_ID)]
  RefExpFeatureMap[, REFERENCE:=as.integer(REFERENCE)]
  data.table::setkeyv(RefExpFeatureMap, c("ALIGNMENT_GROUP_ID", "REFERENCE", "FEATURE_ID"))
  RefExpFeatureMap <- unique(RefExpFeatureMap)

  if (oswMerged){
    oswfile <- unique(fileInfo$featureFile)
    conn <- DBI::dbConnect(RSQLite::SQLite(), oswfile)
    ## Write table to database, overwrite if one already exits
    DBI::dbWriteTable(conn, "ALIGNMENT_GROUP_FEATURE_MAPPING", RefExpFeatureMap, overwrite=TRUE)
    DBI::dbDisconnect(conn)
  } else {
    feature_id_mapping_file <- paste(gsub("[.]tsv", "", outFile), "reference_experiment_feature_map.tsv", sep="_")
    utils::write.table(RefExpFeatureMap, file = feature_id_mapping_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  invisible(NULL)
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

# Alignment of 17186, sequence = 7351_ANAMGIPSLTVTNVPGSTLSR/3 across runs
# ref = hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt
# eXp = hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt
# Example: test_populateReferenceExperimentFeatureAlignmentMap
testtAligned <- function(){
  tAligned <- structure(c(4313, 4316.4, 4319.8, 4323.2, 4326.6, 4330.1, 4333.5,
                          4336.9, 4340.3, 4343.7, 4347.1, 4350.5, 4353.9, 4357.4, 4360.8,
                          4364.2, 4367.6, 4371, 4374.4, 4377.8, 4381.3, 4384.7, 4388.1,
                          4391.5, 4394.9, 4398.3, 4401.7, 4405.2, 4408.6, 4412, 4415.4,
                          4418.8, 4422.2, 4425.6, 4429.1, 4432.5, 4435.9, 4439.3, 4442.7,
                          4446.1, 4449.5, 4452.9, 4456.4, 4459.8, 4463.2, 4466.6, 4470,
                          4473.4, 4476.8, 4480.3, 4483.7, 4487.1, 4490.5, 4493.9, 4497.3,
                          4500.7, 4504.2, 4507.6, 4511, 4514.4, 4517.8, 4521.2, 4524.6,
                          4528.1, 4531.5, 4534.9, 4538.3, 4541.7, 4545.1, 4548.5, 4551.9,
                          4555.4, 4558.8, 4562.2, 4565.6, 4569, 4572.4, 4575.8, 4579.3,
                          4582.7, 4586.1, 4589.5, 4592.9, 4596.3, 4599.7, 4603.2, 4606.6,
                          4610, 4613.4, 4616.8, 4620.2, 4623.6, 4627, 4630.5, 4633.9, 4637.3,
                          4640.7, 4644.1, 4647.5, 4650.9, 4654.4, 4657.8, 4661.2, 4664.6,
                          4668, 4671.4, 4674.8, 4678.3, 4681.7, 4685.1, 4688.5, 4691.9,
                          4695.3, 4698.7, 4702.2, 4705.6, 4709, 4712.4, 4715.8, 4719.2,
                          4722.6, 4726.1, 4729.5, 4732.9, 4736.3, 4739.7, 4743.1, 4746.5,
                          4749.9, 4753.4, 4756.8, 4760.2, 4763.6, 4767, 4770.4, 4773.8,
                          4777.3, 4780.7, 4784.1, 4787.5, 4790.9, 4794.3, 4797.7, 4801.2,
                          4804.6, 4808, 4811.4, 4814.8, 4818.2, 4821.6, 4825.1, 4828.5,
                          4831.9, 4835.3, 4838.7, 4842.1, 4845.5, 4849, 4852.4, 4855.8,
                          4859.2, 4862.6, 4866, 4869.4, 4872.8, 4876.3, 4879.7, 4883.1,
                          4886.5, 4889.9, 4893.3, 4896.7, 4900.2, 4903.6, 4907, NA, 4319.8,
                          4326.6, 4330.1, 4333.5, 4336.9, 4340.3, 4347.1, 4350.5, 4354,
                          4357.4, 4360.8, 4364.2, 4371, 4374.4, 4377.9, 4381.3, 4384.7,
                          4391.5, 4394.9, 4398.3, 4401.8, 4405.2, 4412, 4415.4, 4418.8,
                          4422.2, 4425.6, 4429.1, 4435.9, 4439.3, 4442.7, 4446.1, 4449.5,
                          4456.4, 4459.8, 4463.2, 4466.6, 4470, 4473.4, 4480.3, 4483.7,
                          4487.1, 4490.5, 4493.9, 4500.7, 4504.2, 4507.6, 4511, 4514.4,
                          4517.8, 4524.6, 4528.1, 4531.5, 4534.9, 4538.3, 4545.1, 4548.5,
                          4552, 4555.4, 4558.8, 4562.2, 4565.6, 4569, 4572.4, 4575.9, 4579.3,
                          4582.7, 4586.1, 4589.5, 4592.9, 4596.3, 4599.7, 4603.2, 4604.9,
                          4606.6, 4610, 4613.4, 4616.8, 4620.2, 4623.6, 4627.1, 4630.5,
                          4633.9, 4635.6, 4637.3, 4640.7, 4644.1, 4647.5, 4651, 4654.4,
                          4657.8, 4661.2, 4664.6, 4668, 4671.4, 4673.15, 4674.9, 4678.3,
                          4681.7, 4685.1, 4688.5, 4691.9, 4695.3, 4698.7, 4702.2, 4705.6,
                          4709, 4712.4, 4715.8, 4719.2, 4722.6, 4726.1, 4729.5, 4736.3,
                          4739.7, 4743.1, 4746.5, 4750, 4753.4, 4756.8, 4760.2, 4763.6,
                          4767, 4770.4, 4773.9, 4777.3, 4780.7, 4784.1, 4787.5, 4790.9,
                          4794.3, 4797.8, 4801.2, 4804.6, 4808, 4811.4, 4814.8, 4818.2,
                          4821.6, 4825.1, 4828.5, 4831.9, 4835.3, 4838.7, 4842.1, 4845.5,
                          4849, 4852.4, 4855.8, 4859.2, 4862.6, 4866, 4869.4, 4872.9, 4876.3,
                          4879.7, 4883.1, 4886.5, 4889.9, 4893.3, 4896.8, 4900.2, 4903.6,
                          4907, 4910.4, 4913.8, 4917.2, NA, NA, NA, NA, NA, NA, NA), .Dim = c(175L, 2L))
  tAligned
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

  if(!any(params[["treeDist"]] %in% c("rsquared", "count",  "RSE"))){
    stop("treeDist must be either rsquared, count or RSE.")
  }

  if(!any(params[["treeAgg"]] %in% c("single", "average",  "complete"))){
    stop("treeAgg must be either single, average or complete.")
  }

  if(params[["prefix"]] == "run"){
    stop("prefix cannot be run.")
  }
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

  if(!params[["criterion"]] %in% 1:4){
    stop("criterion must be from 1,2,3,4. Recommended value is 2")
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
#' \item{treeDist}{(string) the method used to build distance matrix. Must be either "rsquared", "count" or "RSE".}
#' \item{treeAgg}{(string) the method used for agglomeration while performing hierarchical clustering. Must be either "single", "average" or "complete".}
#' \item{alignToRoot}{(logical) if TRUE, align leaves to the root in hierarchical clustering, else use already save aligned vectors.}
#' \item{prefix}{(string) name to be used to define merged runs.}
#' \item{context}{(string) used in pyprophet peptide. Must be either "run-specific", "experiment-wide", or "global".}
#' \item{unalignedFDR}{(numeric) must be between 0 and maxFdrQuery. Features below unalignedFDR are
#'  considered for quantification even without the RT alignment.}
#' \item{alignedFDR1}{(numeric) must be between unalignedFDR and alignedFDR2. Features below alignedFDR1 and aligned to the reference are
#'  considered for quantification.}
#' \item{alignedFDR2}{(numeric) must be between alignedFDR1 and maxFdrQuery. Features below alignedFDR2 and within certain distance from the aligned time are
#'  considered for quantification after the alignment.}
#' \item{criterion}{(integer) strategy to select peak if found overlapping peaks. 1:intensity, 2: RT overlap, 3: mscore, 4: edge distance}
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
#' \item{goFactor}{(numeric) penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty. Should be between 10-1000.}
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
                  maxFdrQuery = 0.05, maxIPFFdrQuery = 0.05, maxPeptideFdr = 0.01,
                  analyteFDR = 0.01, treeDist = "count", treeAgg = "single", alignToRoot = FALSE, prefix = "master",
                  context = "global", unalignedFDR = 0.00, alignedFDR1 = 0.05, alignedFDR2 = 0.05, criterion = 2L,
                  level = "Peptide", integrationType = "intensity_sum", baselineType = "base_to_base", fitEMG = FALSE,
                  recalIntensity = FALSE, fillMissing = TRUE, baseSubtraction = FALSE,
                  XICfilter = "sgolay", polyOrd = 4L, kernelLen = 11L,
                  globalAlignment = "loess", globalAlignmentFdr = 0.01, globalAlignmentSpan = 0.1,
                  RSEdistFactor = 3.5, normalization = "mean", simMeasure = "dotProductMasked",
                  alignType = "hybrid", goFactor = 0.125, geFactor = 40,
                  cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                  dotProdThresh = 0.96, gapQuantile = 0.5, kerLen = 9,
                  hardConstrain = FALSE, samples4gradient = 1L,
                  wF = base::min, fillMethod = "spline", splineMethod = "natural", mergeTime = "avg", smoothPeakArea = FALSE,
                  keepFlanks = TRUE, batchSize = 1000L, transitionIntensity = FALSE,
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

# helper function to skip tests if we don't have the 'foo' module
skip_if_no_pyopenms <- function() {
  envName <- "TricEnvr"
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

getBestPkIdx <- function(df, pk, idx, criterion = 2L){
  if(criterion == 1L){ # Check based on intensity
    idx <- idx[which.max(totalIntensity(df, idx))]
    return(idx)
  } else if(criterion == 2L){ # Check based on RT overlap
    idx <- idx[which.max(overlapLen(df, pk, idx))]
    return(idx)
  } else if(criterion == 3L){ # Check based on m-score
    mscores <- .subset2(df, "m_score")[idx]
    i <- which.min(mscores)
    if(sum(mscores[i]>=mscores, na.rm=T) == 1){
      return(idx[i])
    }
    idx <- idx[which(mscores[i]>=mscores)]
  }
  # Check based on edge-distance
  idx <- idx[which.min(pkDist(df, pk, idx))]
  idx
}

pkDist <- function(df, pk, idx){
  left <- pk[1] - .subset2(df, "leftWidth")[idx]
  right <- pk[2] - .subset2(df, "rightWidth")[idx]
  pmin(abs(right) , abs(left))
}

overlapLen <- function(df, pk, idx){
  left <- pmax(pk[1], .subset2(df, "leftWidth")[idx])
  right <- pmin(pk[2], .subset2(df, "rightWidth")[idx])
  right - left
}

totalIntensity <- function(df, idx){
  sapply(.subset2(df, "intensity")[idx], sum)
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
  strategy <- params[["treeDist"]]
  message("Calculating distance matrix using ", strategy)
  if(strategy == "rsquared"){
    distMat <- distMat.rSqrd(features, params, applyFun)
  } else if(strategy == "count"){
    distMat <- distMat.count(features, params, applyFun)
  } else if(strategy == "RSE"){
    distMat <- distMat.RSE(features, params, applyFun)
  } else{
    stop("treeDist must be either rsquared, count or RSE.")
  }
  stats::as.dist(distMat, diag = FALSE, upper = FALSE)
}

distMat.RSE <- function(features, params, applyFun = lapply){
  ##### Get distances among runs based on the Residual Standard Error (RSE). #####
  # Uses params[["globalAlignmentFdr"]] for cut-off because this will be used for global alignment.
  runs <- names(features)
  rseMat <- sapply(runs, function(ref){
    unlist(applyFun(runs, function(eXp){
      fit <- getGlobalAlignment(features, ref, eXp,  params[["globalAlignment"]], params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
      getRSE(fit,  params[["globalAlignment"]])
    }))
  }, USE.NAMES = FALSE)
  # rseMat <- rseMat/max(rseMat) # The denominator should be the length of XICs.
  colnames(rseMat) <- rownames(rseMat) <- runs
  rseMat
}

distMat.rSqrd <- function(features, params, applyFun = lapply){
  # Uses params[["globalAlignmentFdr"]] for cut-off because this will be used for global alignment.
  runs <- names(features)
  distMat <- sapply(runs, function(ref){
    c(applyFun(runs, function(eXp){
      RUNS_RT <- tryCatch(expr = getRTdf(features, ref, eXp, params[["globalAlignmentFdr"]]),
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
  # Use params[["analyteFDR"]] instead of params[["globalAlignmentFdr"]] to provide control about features to choose.
  tmp <- applyFun(features, function(df)
    df[df[["m_score"]] <= params[["analyteFDR"]] & df[["peak_group_rank"]] == 1L, "transition_group_id"][[1]])
  #tmp <- tmp[order(names(tmp), decreasing = FALSE)]
  simMat <- crossprod(table(utils::stack(tmp))) # Number of common peaks in each pair
  a <- sapply(tmp, length)
  m <- matrix(rep(a, each = length(a)) + rep(a, times = length(a)), nrow = length(a))
  distMat <- 1-(2*simMat/m) # 1 - 2*nCommon/nA + nB
  distMat
}


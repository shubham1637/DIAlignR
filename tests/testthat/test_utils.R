context("Utility functions")

test_that("test_getRefRun", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide", maxPeptideFdr = 1.00)
  peptideIDs <- unique(precursors$peptide_id)
  peptidesInfo <- getPeptideScores(fileInfo, peptideIDs, oswMerged = TRUE, runType = "DIA_Proteomics", context = "experiment-wide")
  peptidesInfo <- lapply(peptideIDs, function(pep) dplyr::filter(peptidesInfo, .data$peptide_id == pep))
  names(peptidesInfo) <- as.character(peptideIDs)
  outData <- getRefRun(peptidesInfo)

  expData <- data.table("peptide_id" = c(15L, 105L),
                        "run" = c("run1","run0"), key = "peptide_id")
  expect_identical(dim(outData), c(297L, 2L))
  expect_identical(outData[2:3,], expData)
})

test_that("test_getMultipeptide", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide", maxPeptideFdr = 0.05)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  outData <- getMultipeptide(precursors, features)

  expData <- data.table("transition_group_id" = rep(c(9719L, 9720L),3),
                        "feature_id" = bit64::as.integer64(c(NA_real_, "7742456255764097691", NA_real_, NA_real_, "6462000664077079508", "5135268764240690321", NA_real_, NA_real_, NA_real_, "298844719207353347", NA_real_, NA_real_)),
                        "RT" = c(NA_real_, 2541.83, NA_real_, NA_real_, 2586.12, 2585.61, NA_real_, NA_real_, NA_real_, 2607.05, NA_real_, NA_real_),
                        "intensity" = c(NA_real_, 34.7208, NA_real_, NA_real_, 26.2182, 52.9595, NA_real_, NA_real_, NA_real_, 33.5961, NA_real_, NA_real_),
                        "leftWidth" = c(NA_real_, 2526.555, NA_real_, NA_real_, 2571.738, 2564.094, NA_real_, NA_real_, NA_real_, 2591.431, NA_real_, NA_real_),
                        "rightWidth" = c(NA_real_, 2560.693, NA_real_, NA_real_, 2609.288, 2605.060, NA_real_, NA_real_, NA_real_, 2625.569, NA_real_, NA_real_),
                        "peak_group_rank" = c(NA_integer_, 1L, NA_integer_, NA_integer_, 1L, 1L, NA_integer_, NA_integer_, NA_integer_, 1L, NA_integer_, NA_integer_),
                        "m_score" = c(NA_real_, 5.692077e-05, NA_real_, NA_real_, 1.041916e-03, 5.692077e-05, NA_real_, NA_real_, NA_real_, 2.005418e-04, NA_real_, NA_real_),
                        "run" = rep(c("run0", "run1", "run2"), each = 4),
                        "alignment_rank" = c(rep(NA_integer_,6)),
                        key = "run")
  expect_identical(length(outData), 229L)
  expect_equal(outData[[110]], expData, tolerance = 1e-04)
  expect_equal(outData[["9861"]], expData, tolerance = 1e-03)

  # Test IPF
  dataPath <- system.file("ptms", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_IPF", context = "experiment-wide", maxPeptideFdr = 0.05, useIdentifying = FALSE)
  features <- getFeatures(fileInfo, runType = "DIA_IPF", maxFdrQuery = 0.05, maxIPFFdrQuery = 1)
  outData <- getMultipeptide(precursors, features, runType = "DIA_IPF")

  expData <- data.table("transition_group_id" = c(630L, 630L, 631L, 630L, 631L, 630L, 630L, 630L, 631L, 631L, 631L, 630L, 631L, 630L, 630L, 630L, 630L, 631L, 631L, 631L, 631L, 630L, 631L),
                        "feature_id" = bit64::as.integer64(c(-6191423739689805270, 5547167752232279283, -7097443329173487547, NA_real_, NA_real_, -2549840322423605148, -499421479565831329, -8967853951083755886, 6457111977431447159, 162501776421540861, 4822819366799364222, NA_real_, NA_real_, -8985027232552631382, -3426545446930513222, -1477625472540680028, -4378853637469637664, 7193352087667100241, 4324971034249812598, -8641225003924195901, 166050648022998207, NA_real_, NA_real_)),
                        "RT" = c(1843.36, 1907.49, 1839.59, NA_real_, NA_real_, 1963.89, 1889.3, 1991.81, 1890.85, 1989.07, 1966.6, NA_real_, NA_real_, 1835.29, 1865.66, 1769.11, 1715.97, 1769.75, 1836.3, 1866.79, 1717.79, NA_real_, NA_real_),
                        "intensity" = c(4970, 2833, 11270, NA_real_, NA_real_, 8603, 14265, 2476, 23204, 2833, 19753, NA_real_, NA_real_, 51890, 14820, 115025, 7271, 112106, 57985, 13743, 8220, NA_real_, NA_real_),
                        "leftWidth" = c(1816.43005371094, 1885.5, 1814.77001953125, NA_real_, NA_real_, 1943.56994628906, 1863.59997558594, 1983.55004882812, 1861.93994140625, 1981.89001464844, 1938.27001953125, NA_real_, NA_real_, 1820.07995605469, 1852.7900390625, 1747.38000488281, 1700.13000488281, 1745.71997070312, 1818.42004394531, 1854.77001953125, 1698.46997070312, NA_real_, NA_real_),
                        "rightWidth" = c(1860.05004882812, 1921.83996582031, 1872.93005371094, NA_real_, NA_real_, 1983.55004882812, 1910.84997558594, 2008.98999023438, 1912.82995605469, 2003.69995117188, 1985.52001953125, NA_real_, NA_real_, 1852.7900390625, 1881.86999511719, 1787.36999511719, 1736.47998046875, 1789.33996582031, 1854.77001953125, 1883.84997558594, 1734.81994628906, NA_real_, NA_real_),
                        "peak_group_rank" = c(1, 2, 1, NA_integer_, NA_integer_, 1, 2, 3, 1, 2, 3, NA_integer_, NA_integer_, 1, 2, 3, 4, 1, 2, 3, 4, NA_integer_, NA_integer_),
                        "m_score" = c(0.386984247287001, 0.428802116378765, 0.522428735921873, NA_real_, NA_real_, 0.0427722787632167, 0.391141970210945, 0.552582359159911, 0.269746926304544, 0.49419409884214, 0.18446161248398, NA_real_, NA_real_, 0.00344456667388718, 0.149802921806063, 0.744984159931559, 0.342361325716167, 0.710129418011622, 0.00989114312849029, 0.526409498125047, 0.502204883032099, NA_real_, NA_real_),
                        "run" = c('run0', 'run0', 'run0', 'run0', 'run0', 'run1', 'run1', 'run1', 'run1', 'run1', 'run1', 'run1', 'run1', 'run2', 'run2', 'run2', 'run2', 'run2', 'run2', 'run2', 'run2', 'run2', 'run2'),
                        "alignment_rank" = c(rep(NA_integer_,23)),
                        "ms2_m_score" = c(0.00555125306455705, 0.00738654214619651, 0.0449198676982163, NA_real_, NA_real_, 0.00278421980327629, 0.0042352644409684, 0.0449911911702401, 0.0037987931037582, 0.0042352644409684, 0.00643545542248512, NA_real_, NA_real_, 0.002980631985151, 0.0042352644409684, 0.0042352644409684, 0.00925691124810667, 0.00278421980327629, 0.00278421980327629, 0.002980631985151, 0.0042352644409684, NA_real_, NA_real_),
                        key = "run")
  expect_identical(length(outData), 1L)
  expect_equal(outData[[1]], expData, tolerance = 1e-04)
})

test_that("test_writeTables", {
  # To check if the two files are same, we can use tools::md5sum(). However, this will
  # not tell us which value is different.
  dataPath <- system.file("extdata", package = "DIAlignR")
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  # data(multipeptide_DIAlignR, package="DIAlignR")
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide")
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  allIDs <- unique(unlist(lapply(features, function(df) df[m_score <= 0.05, transition_group_id]),
                          recursive = FALSE, use.names = FALSE))
  precursors <- precursors[data.table::CJ(unique(peptide_id), allIDs) , nomatch=0L]
  multipeptide <- getMultipeptide(precursors, features)
  outData <- writeTables(fileInfo, multipeptide, precursors)
  expect_identical(nrow(outData), 500L)

  multipeptide[["7040"]][c(1,3), alignment_rank:= 1L]
  multipeptide[["3200"]][c(1,3, 5), alignment_rank:= 1L]
  outData <- writeTables(fileInfo, multipeptide, precursors)

  expect_identical(dim(outData), c(502L, 14L))
  expData <- data.table(peptide_id = c(3200L, 7040L), precursor = c(523L, 32L),
                        run = fileInfo$runName[1:2], RT = c(1529.54, NA_real_),
                        intensity = c(13.5686, NA_real_), leftWidth = c(1512.653,NA_real_), rightWidth= c(1546.791,NA_real_),
                        peak_group_rank = c(1L, NA_integer_), m_score = c(0.01440752, NA_real_), alignment_rank = 1L,
                        feature_id = c("2382660248384924660", NA_character_), sequence = c("DVDQYPR", "GNNSVYMNNFLNLILQNER"),
                        charge = c(2L, 3L), group_label = c("10483_DVDQYPR/2", "10030_GNNSVYMNNFLNLILQNER/3"),
                        key = NULL)
  expect_equal(outData[c(82,167),], expData, tolerance = 1e-05)
})

test_that("test_checkParams", {
})

test_that("test_paramsDIAlignR", {
})

test_that("test_alignmentStats", {
})


test_that("test_checkOverlap",{
  expect_true(checkOverlap(c(1.1, 3.1), c(2.1, 3.1)))
  expect_true(checkOverlap(c(1.1, 3.1), c(3.1, 4.1)))
  expect_true(checkOverlap(c(1.1, 3.1), c(2.1, 2.5)))
  expect_true(checkOverlap(c(1.1, 3.1), c(0.1, 9.1)))
  expect_false(checkOverlap(c(1.1, 3.1), c(3.2, 7.1)))

})

test_that("test_ipfReassignFDR", {
  dataPath <- system.file("ptms", package = "DIAlignR")
  params <- paramsDIAlignR()
  params$chromFile <- "sqMass"
  params$runType <- "DIA_IPF"
  params$maxFdrQuery <- 0.05
  params$maxIPFFdrQuery <- 1
  params$useIdentifying <- TRUE
  params$unalignedFDR <- 0.05
  params$alignedFDR2 <- 0.05
  params$globalAlignment <- "linear"
  params$samples4gradient <- 100L
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
  precursors <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_IPF", context = "experiment-wide", maxPeptideFdr = 0.05, useIdentifying = FALSE)
  peptideIDs <- precursors[, logical(1), keyby = peptide_id]$peptide_id
  peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged=TRUE, runType = "DIA_IPF", context ="global")
  refRuns <- data.table(peptide_id=174L, run="run2")
  features <- getFeatures(fileInfo, runType = "DIA_IPF", maxFdrQuery = 0.05, maxIPFFdrQuery = 1)
  mzPntrs <- getMZMLpointers(fileInfo)
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs, lapply)
  multipeptide <- getMultipeptide(precursors, features, runType = "DIA_IPF")
  globalFits <- getGlobalFits(refRuns, features, fileInfo, "linear", 0.01, 0.1, lapply)
  RSE <- lapply(globalFits, getRSE, "linear")
  globalFits <- lapply(globalFits, extractFit, "linear")
  perBatch(iBatch=1, peptideIDs, multipeptide, refRuns, precursors, prec2chromIndex, fileInfo, mzPntrs, params, globalFits, RSE, applyFun = lapply)
  ## Clean up
  for(mz in mzPntrs){
    if(is(mz)[1] == "SQLiteConnection") DBI::dbDisconnect(mz)
    if(is(mz)[1] == "mzRpwiz") rm(mz)
  }
  rm(prec2chromIndex, globalFits, RSE)
  rm(features)

  finalTbl <- writeTables(fileInfo, multipeptide, precursors)
  finalTbl <- ipfReassignFDR(dt = finalTbl, refRuns, fileInfo, params)

  expData <- data.table( peptide_id = rep(c(174L), 6),
                         precursor = c(630L, 630L, 630L, 631L, 631L, 631L),
                         run = c('chludwig_K150309_004b_SW_1_16', 'chludwig_K150309_008_SW_1_4', 'chludwig_K150309_013_SW_0', 'chludwig_K150309_004b_SW_1_16', 'chludwig_K150309_008_SW_1_4', 'chludwig_K150309_013_SW_0'),
                         RT = c(1909.1, 1963.89, 1835.29, 1909.1, 1966.6, 1836.3),
                         intensity = c(1142.37258827012, 8603, 51890, 4651.38812575098, 19753, 57985),
                         leftWidth = c(1894.6, 1943.56994628906, 1820.07995605469, 1894.6, 1938.27001953125, 1818.42004394531),
                         rightWidth = c(1927.3, 1983.55004882812, 1852.7900390625, 1927.3, 1985.52001953125, 1854.77001953125),
                         peak_group_rank = c(NA_integer_, 1L, 1L, NA_integer_, 3L, 2L),
                         original_m_score = c(NA_real_, 0.0427722787632167, 0.00344456667388718, NA_real_, 0.18446161248398, 0.00989114312849029),
                         alignment_rank = c(1L, 1L, 1L, 1L, 1L, 1L),
                         feature_id = c(NA_character_, "-2549840322423605148", "-8985027232552631382", NA_character_, "4822819366799364222", "4324971034249812598"),
                         sequence = c('ANS(UniMod:21)SPTTNIDHLK(UniMod:259)', 'ANS(UniMod:21)SPTTNIDHLK(UniMod:259)', 'ANS(UniMod:21)SPTTNIDHLK(UniMod:259)', 'ANS(UniMod:21)SPTTNIDHLK(UniMod:259)', 'ANS(UniMod:21)SPTTNIDHLK(UniMod:259)', 'ANS(UniMod:21)SPTTNIDHLK(UniMod:259)'),
                         charge = c(2L, 2L, 2L, 3L, 3L, 3L),
                         group_label = c('1565', '1565', '1565', '1566', '1566', '1566'),
                         ms2_m_score = c(NA_real_, 0.00278421980327629, 0.002980631985151, NA_real_, 0.00643545542248512, 0.00278421980327629),
                         ref_run = c('chludwig_K150309_013_SW_0', 'chludwig_K150309_013_SW_0', 'chludwig_K150309_013_SW_0', 'chludwig_K150309_013_SW_0', 'chludwig_K150309_013_SW_0', 'chludwig_K150309_013_SW_0'),
                         m_score = c(0.05, 0.00278421980327629, 0.00344456667388718, 0.05, 0.00643545542248512, 0.00989114312849029)
                       )
  setkey(expData, "peptide_id")

  expect_identical(dim(finalTbl), c(6L, 17L))
  expect_equal(finalTbl, expData, tolerance = 1e-05)
})

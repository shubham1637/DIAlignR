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

test_that("test_getRefExpFeatureMap",
          {
            dataPath <- system.file("extdata", package = "DIAlignR")
            fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
            precursors <- getPrecursors(fileInfo, oswMerged = TRUE, context = "experiment-wide", maxPeptideFdr = 0.05)
            features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
            outData <- getRefExpFeatureMap(precursors, features)
            expData <- data.table("reference_feature_id" = rep(0L, 15), "experiment_feature_id" = rep(0L, 15))
            setattr(expData[['reference_feature_id']], "class","integer64")
            setattr(expData[['experiment_feature_id']], "class","integer64")
            setkeyv(expData, "reference_feature_id")
            expect_identical(length(outData), 229L)
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

test_that("test_writeOutFeatureAlignmentMap",
          {
            #### Prepare data and outout ####
            dataPath <- system.file("extdata", package = "DIAlignR")
            fileInfo <- getRunNames(dataPath, oswMerged = TRUE)

            multiFeatureAlignmentMap <- list(`1338` = structure(list(reference_feature_id = structure(c(4.75292059933946e+94,
                                                                                                        5.52612460272647e-11, 3.59272034753732e+234, 3.74370892823184e+230,
                                                                                                        2.34380646139584e+234, 4.75292059933946e+94, 5.52612460272647e-11,
                                                                                                        3.59272034753732e+234, 3.74370892823184e+230, 0, 0, 0, 0, 0,
                                                                                                        0), class = "integer64"),
                                                                     experiment_feature_id = structure(c(2.04161779709829e+87,
                                                                                                         1.60545397951231e+66, 5.20548354836881e-89, 4.57942078616298e+306,
                                                                                                         2.66188671607728e-167, 1.89738801704963e+63, 7.46666102909777e-250,
                                                                                                         1.01955111889378e-75, 5.12280305276836e+252, 0, 0, 0, 0, 0, 0
                                                                                                        ), class = "integer64")), row.names = c(NA, -15L), class = c("data.table", "data.frame")),
                                             `19045` = structure(list(reference_feature_id = structure(c(6.28056556040491e-257,
                                                                                                         4.37872094299259e+26, 4.97487068543727e-50, 1.96129357589358e+225,
                                                                                                         6.28056556040491e-257, 4.37872094299259e+26, 4.97487068543727e-50,
                                                                                                         2.22889297151144e-67, 1.96129357589358e+225, 0, 0, 0, 0,
                                                                                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), class = "integer64"),
                                                                      experiment_feature_id = structure(c(1.39189055078857e+239,
                                                                                                          1.64937896938207e+218, 1.3591112627868e-271, 3.36180583288588e-229,
                                                                                                          2.16310643042558e-237, 5.82423961711318e+254, 9.33806679439292e+229,
                                                                                                          1.58730393504064e-226, 7.21410639770476e+159, 0, 0, 0,
                                                                                                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), class = "integer64")), row.names = c(NA, -30L), class = c("data.table", "data.frame")),
                                             `19800` = structure(list(reference_feature_id = structure(c(4.72799309850901e-140,
                                                                                                         1.98100991301857e-234, 4.7297745210988e+87, 4.72799309850901e-140,
                                                                                                         4.7297745210988e+87, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), class = "integer64"),
                                                                      experiment_feature_id = structure(c(8.56560434196891e+243,
                                                                                                          1.62593959656114e+49, 4.1149805476014e-97, 1.94083176316202e-40,
                                                                                                          1.29045098420263e-28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), class = "integer64")), row.names = c(NA, -15L), class = c("data.table", "data.frame")))

            #### Expected Data ####
            expData <- structure(list(ALIGNMENT_GROUP_ID = c(1L, 1L, 1L, 2L, 2L, 3L,
                                                             3L, 3L, 4L, 4L, 5L, 5L, 5L, 6L, 6L, 6L, 7L, 7L, 7L, 8L, 8L, 8L,
                                                             9L, 9L, 9L, 10L, 10L, 10L, 11L, 11L, 11L, 12L, 12L, 13L, 13L,
                                                             13L),
                                      REFERENCE = c(0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 0L,
                                                    0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L,
                                                    1L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L),
                                      FEATURE_ID = structure(c(2.16310643042558e-237,
                                                               1.39189055078857e+239, 6.28056556040491e-257, 1.62593959656114e+49,
                                                               1.98100991301857e-234, 1.94083176316202e-40, 8.56560434196891e+243,
                                                               4.72799309850901e-140, 1.58730393504064e-226, 2.22889297151144e-67,
                                                               1.3591112627868e-271, 9.33806679439292e+229, 4.97487068543727e-50,
                                                               7.46666102909777e-250, 1.60545397951231e+66, 5.52612460272647e-11,
                                                               1.64937896938207e+218, 5.82423961711318e+254, 4.37872094299259e+26,
                                                               4.1149805476014e-97, 1.29045098420263e-28, 4.7297745210988e+87,
                                                               1.89738801704963e+63, 2.04161779709829e+87, 4.75292059933946e+94,
                                                               3.36180583288588e-229, 7.21410639770476e+159, 1.96129357589358e+225,
                                                               5.12280305276836e+252, 4.57942078616298e+306, 3.74370892823184e+230,
                                                               2.66188671607728e-167, 2.34380646139584e+234, 5.20548354836881e-89,
                                                               1.01955111889378e-75, 3.59272034753732e+234), class = "integer64")),
                                 class = "data.frame", row.names = c(NA, 36L))

            #### Writing to OSW ####
            # Make a Copy of fileInfo with a copy of the OSW file to avoid git history of test merged.osw file from being changed
            fileInfotmp <- fileInfo
            copiedOSWpath <- paste(getwd(), basename(unique(fileInfo$featureFile)[1]), sep=.Platform$file.sep)
            file.copy(unique(fileInfo$featureFile)[1], copiedOSWpath)
            fileInfotmp$featureFile <- copiedOSWpath

            # Write out alignment map to disk: OSW
            writeOutFeatureAlignmentMap(multiFeatureAlignmentMap, oswMerged=TRUE, fileInfotmp)

            # Load written out table from osw file
            con <- DBI::dbConnect(RSQLite::SQLite(), fileInfotmp$featureFile[1])
            outData <- DBI::dbGetQuery(con, "SELECT * FROM ALIGNMENT_GROUP_FEATURE_MAPPING")
            DBI::dbDisconnect(con)

            #### Test Expectations ####
            expect_equal(outData, expData, tolerance = 1e-6)

            #### Clean up and Remove Test Output ####
            file.remove(copiedOSWpath)

            #### Writing to TSV ####
            # Write out alignment map to disk: TSV
            writeOutFeatureAlignmentMap(multiFeatureAlignmentMap, oswMerged=FALSE, fileInfotmp)

            # Load written out table from TSV file
            outData <- read.table("reference_experiment_feature_map.tsv", header = TRUE, sep = "\t", colClasses = c("integer", "integer", "integer64"))

            #### Test Expectations ####
            expect_equal(outData, expData, tolerance = 1e-6)

            #### Clean up and Remove Test Output ####
            file.remove("reference_experiment_feature_map.tsv")
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
  params$alignedFDR1 <- 0.05
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
  RSE <- lapply(globalFits, getRSE, params$globalAlignment)
  globalFits <- lapply(globalFits, extractFit, params$globalAlignment)
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
                         RT = c(NA_real_, 1963.89, 1835.29, NA_real_, 1966.6, 1836.3),
                         intensity = c(NA_real_, 8603, 51890, NA_real_, 19753, 57985),
                         leftWidth = c(NA_real_, 1943.56994628906, 1820.07995605469, NA_real_, 1938.27001953125, 1818.42004394531),
                         rightWidth = c(NA_real_, 1983.55004882812, 1852.7900390625, NA_real_, 1985.52001953125, 1854.77001953125),
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

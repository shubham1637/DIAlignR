context("Minimum Spanning Tree")

test_that("test_getMST",{
  m <- matrix(c(0,1,2,3, 1,0,1.5,1.5, 2,1.5,0,1, 3,1.5,1,0), byrow = TRUE,
               ncol = 4, dimnames = list(c("run1", "run2", "run3", "run4"),
                                       c("run1", "run2", "run3", "run4")))
  distMat <- as.dist(m)
  outData <- getMST(distMat)
  expect_identical(outData, cbind("A"=c("run1","run2","run3"),
                              "B"=c("run2","run3","run4")))
})

test_that("test_traverseMST",{
  net <- cbind("A"=c("run1","run2","run3"),"B"=c("run2","run3","run4"))
  outData <- traverseMST(net, "run2")
  expect_identical(outData, cbind("ref"=c("run2","run2","run3"),
                                  "eXp"=c("run1","run3","run4")))
  outData <- traverseMST(net, "run4")
  expect_identical(outData, cbind("ref"=c("run4","run3","run2"),
                                  "eXp"=c("run3","run2","run1")))
})

test_that("test_mstAlignRuns",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["baseSubtraction"]] <- TRUE
  params[["context"]] <- "experiment-wide"
  peps <- c(8496L, 15879L, 19800L)
  params[["transitionIntensity"]] <- TRUE
  tree <- "run2 run2\nrun1 run0"
  mstAlignRuns(dataPath = dataPath,  outFile = "temp", params = params, oswMerged = TRUE,
               runs = NULL, peps = peps, mstNet = tree, applyFun = lapply)
  outData <- read.table("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expData <- read.table("test4.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  expect_equal(outData[["RT"]], expData[["RT"]], tolerance = 1e-04)
  x <- sapply(outData[["intensity"]], function(a) as.numeric(strsplit(a, split = ",")[[1]]), USE.NAMES = FALSE)
  y <- sapply(expData[["intensity"]], function(a) as.numeric(strsplit(a, split = ",")[[1]]), USE.NAMES = FALSE)
  expect_equal(x, y, tolerance = 1e-04)
  for(i in 6:10){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.tsv")

  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["baseSubtraction"]] <- TRUE
  params[["maxPeptideFdr"]] <- 0.05
  params[["XICfilter"]] <- "none"
  params[["globalAlignment"]] <- "loess"
  params[["context"]] <- "experiment-wide"
  params[["treeDist"]] <- "rsquared"
  expect_message(
    mstAlignRuns(dataPath = dataPath,  outFile = "temp", params = params, oswMerged = TRUE,
                      runs = NULL, applyFun = lapply)
  )
  outData <- read.table("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expData <- read.table("test5.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 4:10){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.tsv")
})

test_that("test_MSTperBatch",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["unalignedFDR"]] <- 0.01
  params[["baseSubtraction"]] <- TRUE
  params[["context"]] <- "experiment-wide"
  params[["kernelLen"]] <- 13L
  params[["globalAlignmentFdr"]] <- 0.05
  params[["maxPeptideFdr"]] <- 0.05
  params$batchSize <- 1L
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE, params)
  precursors <- getPrecursors(fileInfo, oswMerged= TRUE, params[["runType"]], params[["context"]], params[["maxPeptideFdr"]])
  precursors <- precursors[precursors$peptide_id %in% c("7040", "9861", "14383"),]
  peptideIDs <-  c(7040L, 9861L, 14383L)
  mzPntrs <- getMZMLpointers(fileInfo)
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  multipeptide <- getMultipeptide(precursors, features)
  refRuns <- data.frame("peptide_id" = c("7040", "9861", "14383"), "run" = "run1")
  globalFits <- getGlobalFits(refRuns, features, fileInfo, params[["globalAlignment"]],
                              params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
  RSE <- list()
  RSE[["run1_run2"]] <- RSE[["run1_run0"]] <- 38.6594179136227/params$RSEdistFactor
  globalFits <- lapply(globalFits, extractFit, params[["globalAlignment"]])

  nets <- list(); nets[["run1"]] <- cbind(ref = c("run1", "run2"), eXp = c("run2", "run0"))
  # Case 1 Missing chromatograms
  df <- data.table::copy(multipeptide[["7040"]])
  expect_message(MSTperBatch(iBatch = 1L, nets, peptideIDs, multipeptide, refRuns, precursors,
                          prec2chromIndex, fileInfo, mzPntrs, params, globalFits, RSE))
  expect_equal(df, multipeptide[["7040"]])

  # Case 2 Features are already present. No alignment needed.
  df <- data.table::copy(multipeptide[["14383"]])
  MSTperBatch(iBatch = 3L, nets, peptideIDs, multipeptide, refRuns, precursors, prec2chromIndex, fileInfo,
           mzPntrs, params, globalFits, RSE)
  df$alignment_rank[c(1,3,5)] <- 1L
  expect_equal(multipeptide[["14383"]], df)

  # Case 3
  df <- data.table::copy(multipeptide[["9861"]])
  MSTperBatch(iBatch = 2L, nets, peptideIDs, multipeptide, refRuns, precursors,
           prec2chromIndex, fileInfo, mzPntrs, params, globalFits, RSE)

  data.table::set(df, 1L, c(3L,4L,5L,6L), list(2541.83, 12.92301, 2526.555, 2560.693))
  data.table::set(df, 9L, c(3L,4L,5L,6L), list(2607.05, 11.80541, 2591.431, 2625.569))
  df$alignment_rank[c(1,2,5,6,9,10)] <- 1L
  expect_equal(df, multipeptide[["9861"]], tolerance = 1e-06)
  for(con in mzPntrs) DBI::dbDisconnect(con)
})

test_that("test_alignToRefMST",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["unalignedFDR"]] <- 0.01
  params[["baseSubtraction"]] <- TRUE
  params[["maxPeptideFdr"]] <- 0.05
  params[["context"]] <- "experiment-wide"
  params$kernelLen <- 13L
  params[["globalAlignment"]] <- "linear"
  params[["globalAlignmentFdr"]] <- 0.05
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE, params)
  precursors <- getPrecursors(fileInfo, oswMerged= TRUE, params[["runType"]], params[["context"]], params[["maxPeptideFdr"]])
  precursors <- precursors[precursors$peptide_id %in% c("7040", "9861", "14383"),]

  mzPntrs <- getMZMLpointers(fileInfo)
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  refRuns <- data.table("peptide_id" = c("7040", "14383", "9861"), "run" = "run1", key = "peptide_id")
  net <- cbind(ref = c("run1", "run2"), eXp = c("run2", "run0"))
  globalFits <- getGlobalFits(refRuns, features, fileInfo, params[["globalAlignment"]],
                              params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
  RSE <- list()
  RSE[["run1_run2"]] <- 38.6594179136227/params$RSEdistFactor
  globalFits <- lapply(globalFits, extractFit, params[["globalAlignment"]])

  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  data(multipeptide_DIAlignR, package="DIAlignR")

  # Case 1
  df <- data.table::copy(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[3] <- 1L; df$m_score[5] <- 0.06
  XICs <- list(); XICs[["run2"]] <- list(); XICs[["run1"]] <- list()
  XICs[["run1"]][["4618"]] <- lapply(XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]], as.matrix)
  XICs[["run2"]][["4618"]] <- lapply(XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]], as.matrix)
  alignToRefMST(iNet = 1L, net = net, fileInfo, XICs, params, analytes = 4618L, df, globalFits, RSE)
  expect_equal(df[6,], data.table("transition_group_id" = 4618L, feature_id = bit64::NA_integer64_,
                                  RT = 5241.30, intensity = 189.304,  leftWidth = 5224.2, rightWidth = 5265.2, peak_group_rank = NA_integer_,
                                  m_score = NA_real_, run = "run2", alignment_rank = 1, key = "run"),
               tolerance = 1e-06)

  # Case 2
  df <- data.table::copy(multipeptide_DIAlignR[["14383"]])
  df$alignment_rank[3] <- 1L
  alignToRefMST(iNet = 1L, net = net, fileInfo, NULL, params, analytes = 4618L, df, globalFits, RSE)
  expect_equal(df[c(5,6), alignment_rank], c(1L, NA_integer_))

  # Case 3
  params[["chromFile"]] <- "mzML"
  chromIndices <- prec2chromIndex[["run1"]][c(2,3), chromatogramIndex]
  mz <- mzR::openMSfile(file.path(dataPath, "xics","hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"))
  XICs.ref <- lapply(chromIndices, function(i) extractXIC_group(mz, chromIndices = i+1)) # sqMass file indices start with 0, mzML start with 1
  names(XICs.ref) <- c("9719", "9720")
  chromIndices <- prec2chromIndex[["run2"]][c(2,3), chromatogramIndex]
  mz <- mzR::openMSfile(file.path(dataPath, "xics","hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"))
  xics <- lapply(chromIndices, function(i) extractXIC_group(mz, chromIndices = i+1)) # sqMass file indices start with 0, mzML start with 1
  names(xics) <- c("9719", "9720")
  XICs <- list(); XICs[["run2"]] <- xics; XICs[["run1"]] <- XICs.ref
  df <- multipeptide_DIAlignR[["9861"]]
  df$alignment_rank[6] <- 1L
  alignToRefMST(iNet = 1L, net = net, fileInfo, XICs, params, analytes = c(9719, 9720), df, globalFits, RSE)

  expect_equal(df[10,alignment_rank],1L)
  expect_equal(df[9,], data.table("transition_group_id" = 9719L, feature_id = bit64::NA_integer64_,
                                  RT = 2607.05, intensity = 11.80541,  leftWidth = 2591.431, rightWidth = 2625.569, peak_group_rank = NA_integer_,
                                  m_score = NA_real_, run = "run2", alignment_rank = 1, key = "run"),
               tolerance = 1e-06)
  for(con in mzPntrs) DBI::dbDisconnect(con)
})

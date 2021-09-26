context("Progressive alignment")

test_that("test_progAlignRuns", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["kernelLen"]] <- 9L
  params[["maxPeptideFdr"]] <- 0.05
  params[["context"]] <- "experiment-wide"
  params[["globalAlignment"]] <- "linear"
  params$treeDist <- "rsquared"
  params[["baseSubtraction"]] <- TRUE
  params[["unalignedFDR"]] <- 0.01
  BiocParallel::register(BiocParallel::MulticoreParam())
  params[["chromFile"]] <- "sqMass"
  progAlignRuns(dataPath, params = params, outFile = "temp", ropenms = NULL, applyFun = BiocParallel::bplapply)
  outData <- data.table::fread("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE, integer64 = "integer64")
  expData <- data.table::fread("test6.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE, integer64 = "integer64")
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 4:14) expect_equal(outData[[i]], expData[[i]], tolerance = 1e-05)
  file.remove("temp.tsv")
  file.remove(file.path(dataPath, "master.merged.osw"))
  file.remove(file.path(dataPath, "temp.temp.RData"))
  file.remove(list.files(dataPath, pattern = "*_av.rds", full.names = TRUE))
  file.remove(list.files(file.path(dataPath, "xics"), pattern = "^master[0-9]+\\.chrom\\.sqMass$", full.names = TRUE))
})

test_that("test_progAlignRuns2", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["baseSubtraction"]] <- TRUE
  params[["unalignedFDR"]] <- 0.01
  params[["kernelLen"]] <- 9L
  params[["maxPeptideFdr"]] <- 0.05
  params[["context"]] <- "experiment-wide"
  params[["globalAlignment"]] <- "linear"
  params[["chromFile"]] <- "sqMass"
  params[["transitionIntensity"]] <- TRUE
  params[["fractionNum"]] <- 2L
  params[["fraction"]] <- 1L
  text1 <- "(run1:0.08857142857,(run0:0.06857142857,run2:0.06857142857)masterB:0.02)master1;"
  progTree1(dataPath, params, outFile = "temp", newickTree = text1)
  expect_warning(progSplit2(dataPath, params, outFile = "temp"))
  params[["fraction"]] <- 2L
  progSplit2(dataPath, params, outFile = "temp")
  expect_warning(progComb3(dataPath, params, outFile = "temp"))
  params[["fraction"]] <- 1L
  progSplit4(dataPath, params, outFile = "temp")
  params[["fraction"]] <- 2L
  progSplit4(dataPath, params, outFile = "temp")

  file.remove(list.files(dataPath, pattern = "temp.*rds", full.names = TRUE))
  file.remove(list.files(file.path(dataPath, "xics"), pattern = "^master[A-Za-z0-9]+\\.chrom\\.sqMass$", full.names = TRUE))
  file.remove(list.files(dataPath, pattern = "temp.*RData", full.names = TRUE))
  file.remove(list.files(dataPath, pattern = "*_av.rds", full.names = TRUE))

  df1 <- data.table::fread("temp_1_2.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE, integer64 = "integer64")
  df2 <- data.table::fread("temp_2_2.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE, integer64 = "integer64")
  outData <- data.table::rbindlist(list(df1, df2))
  data.table::setorder(outData, peptide_id, precursor, run)
  outData$intensity <- sapply(outData$intensity, function(x){
    suppressWarnings(y <- as.numeric(strsplit(x, split = ",")[[1]]))
    y[is.nan(y)] <- NA_real_
    sum(y)})
  expData <- data.table::fread("test3.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE, integer64 = "integer64")
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 4:14) expect_equal(outData[[i]], expData[[i]], tolerance = 1e-05)
  file.remove("temp_1_2.tsv")
  file.remove("temp_2_2.tsv")
})

test_that("test_alignToRoot4", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["kernelLen"]] <- 9L
  params[["alignToRoot"]] <- TRUE
  params[["context"]] <- "experiment-wide"
  params[["globalAlignment"]] <- "linear"
  params[["transitionIntensity"]] <- TRUE
  params[["fractionNum"]] <- 2L
  params[["fraction"]] <- 1L
  categ <- data.frame(run = paste0("run", 0:2), ca= c("a", "b", "a"))
  text1 <- "(run1:0.08857142857,(run0:0.06857142857,run2:0.06857142857)master2:0.02)master1;"
  progTree1(dataPath, params, outFile = "temp", categ =  categ, newickTree = text1)
  expect_warning(progSplit2(dataPath, params, outFile = "temp"))
  params[["fraction"]] <- 2L
  progSplit2(dataPath, params, outFile = "temp")
  expect_warning(progComb3(dataPath, params, outFile = "temp"))
  params[["fraction"]] <- 1L
  alignToRoot4(dataPath, params, outFile = "temp")
  params[["fraction"]] <- 2L
  alignToRoot4(dataPath, params, outFile = "temp")
  file.remove(list.files(dataPath, pattern = "temp.*rds", full.names = TRUE))
  file.remove(list.files(file.path(dataPath, "xics"), pattern = "^master[A-Za-z0-9]+\\.chrom\\.sqMass$", full.names = TRUE))
  file.remove(list.files(dataPath, pattern = "temp.*RData", full.names = TRUE))
  file.remove(list.files(dataPath, pattern = "*_av.rds", full.names = TRUE))

  df1 <- data.table::fread("temp_1_2.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE, integer64 = "integer64")
  df2 <- data.table::fread("temp_2_2.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE, integer64 = "integer64")
  outData <- data.table::rbindlist(list(df1, df2))
  data.table::setorder(outData, peptide_id, precursor, run)
  expData <- data.table::fread("test7.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE, integer64 = "integer64")
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  expect_identical(outData[["intensity"]], expData[["intensity"]])
  for(i in c(4,6:14)) expect_equal(outData[[i]], expData[[i]], tolerance = 1e-05)
  file.remove("temp_1_2.tsv")
  file.remove("temp_2_2.tsv")
})

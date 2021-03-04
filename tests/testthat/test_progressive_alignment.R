context("Progressive alignment")

test_that("test_progAlignRuns", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["kernelLen"]] <- 9L
  params[["maxPeptideFdr"]] <- 0.05
  params[["context"]] <- "experiment-wide"
  params[["globalAlignment"]] <- "linear"
  BiocParallel::register(BiocParallel::MulticoreParam())
  params[["chromFile"]] <- "sqMass"
  progAlignRuns(dataPath, params = params, outFile = "temp", ropenms = NULL, applyFun = BiocParallel::bplapply)
  outData <- data.table::fread("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE, integer64 = "integer64")
  expData <- data.table::fread("test3.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE, integer64 = "integer64")
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
  skip_if_no_pyopenms()
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["kernelLen"]] <- 9L
  ropenms <- get_ropenms(condaEnv =  envName)
  params[["context"]] <- "experiment-wide"
  params[["globalAlignment"]] <- "linear"
  BiocParallel::register(BiocParallel::MulticoreParam())
  params[["chromFile"]] <- "mzML"
  for(fun in c(lapply)){
    expect_warning(progAlignRuns(dataPath, params = params, outFile = "temp", ropenms = ropenms, applyFun = fun))
    outData <- data.table::fread("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
    expData <- data.table::fread("test3.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
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
    file.remove(list.files(file.path(dataPath, "xics"), pattern = "^master[0-9]+\\.chrom\\.mzML$", full.names = TRUE))
  }
})

test_that("test_progAlignRuns2", {
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["kernelLen"]] <- 9L
  params[["maxPeptideFdr"]] <- 0.05
  params[["context"]] <- "experiment-wide"
  params[["globalAlignment"]] <- "linear"
  params[["chromFile"]] <- "sqMass"
  params[["transitionIntensity"]] <- TRUE
  params[["fractionNum"]] <- 2L
  params[["fraction"]] <- 1L
  progTree1(dataPath, params, outFile = "temp")
  expect_warning(progSplit2(dataPath, params, outFile = "temp"))
  params[["fraction"]] <- 2L
  progSplit2(dataPath, params, outFile = "temp")
  expect_warning(progComb3(dataPath, params, outFile = "temp"))
  params[["fraction"]] <- 1L
  progSplit4(dataPath, params, outFile = "temp")
  params[["fraction"]] <- 2L
  progSplit4(dataPath, params, outFile = "temp")

  file.remove(list.files(dataPath, pattern = "temp.*rds", full.names = TRUE))
  file.remove(list.files(file.path(dataPath, "xics"), pattern = "^master[0-9]+\\.chrom\\.sqMass$", full.names = TRUE))
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


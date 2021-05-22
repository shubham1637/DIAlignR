# ids <- as.integer(scan(file = "data/ids.txt")) # Peptide Ids manually annotated
# params <- paramsDIAlignR()
# params[["maxFdrQuery"]] <- 1.0
# params[["maxPeptideFdr"]] <- 1.0
# params[["kernelLen"]] <- 11L
# params$globalAlignmentFdr <- 0.001
# params[["batchSize"]] <- 1000L
# BiocParallel::register(BiocParallel::MulticoreParam(workers = 6, log = FALSE, threshold = "INFO", stop.on.error = TRUE))
# alignTargetedRuns(dataPath = ".", params = params, outFile = "star", peps = ids, applyFun = BiocParallel::bplapply)
# mstAlignRuns(dataPath = ".", params = params, outFile = "mst", peps = ids, applyFun = BiocParallel::bplapply)
# ids <- union(ids, as.integer(scan(file = "data/prog.ids.txt")))
# progAlignRuns(dataPath = ".", params = params, outFile = "prog", ropenms = NULL, peps = ids, applyFun = BiocParallel::bplapply)

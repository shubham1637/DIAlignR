#' Peptide quantification through progressive alignment
#'
#' This function expects osw and xics directories at dataPath. It first reads osw files and fetches
#'  chromatogram indices for each analyte. To perform alignment, first a crude guide-tree is built which
#'  can also be provided with newickTree parameter. As we traverse from the leaf-nodes to the root node,
#'  runs are aligned pairwise. The root node is named master1 that has average of all fragment ion chromatograms
#'  and identified peak-groups. These features are propagated back to leaf nodes and finally aligned
#'  features are written in the output file.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-10
#' @importFrom data.table data.table setkeyv rbindlist
#' @inheritParams alignTargetedRuns
#' @param ropenms (pyopenms module) get this python module through \code{\link{get_ropenms}}. Required only for chrom.mzML files.
#' @param newickTree (string) guidance tree in newick format. Look up \code{\link{getTree}}.
#' @return (None)
#' @seealso \code{\link{alignTargetedRuns}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' params[["context"]] <- "experiment-wide"
#' \dontrun{
#' ropenms <- get_ropenms(condaEnv = "envName")
#' progAlignRuns(dataPath, params = params, outFile = "test3", ropenms = ropenms)
#' # Removing aligned vectors
#' file.remove(list.files(dataPath, pattern = "*_av.rds", full.names = TRUE))
#' # Removing temporarily created master chromatograms
#' file.remove(list.files(file.path(dataPath, "xics"), pattern = "^master[0-9]+\\.chrom\\.mzML$", full.names = TRUE))
#' file.remove(file.path(dataPath, "test3.temp.RData"))
#' file.remove(file.path(dataPath, "master.merged.osw"))
#' }
#' @export
progAlignRuns <- function(dataPath, params, outFile = "DIAlignR", ropenms = NULL, oswMerged = TRUE,
                          runs = NULL, peps = NULL, newickTree = NULL, applyFun = lapply){
  #### Check if all parameters make sense.  #########
  if(params[["chromFile"]] == "mzML"){
    if(is.null(ropenms)) stop("ropenms is required to write chrom.mzML files.")
  }
  params <- checkParams(params)

  #### Get filenames from .osw file and check consistency between osw and mzML files. #################
  fileInfo <- getRunNames(dataPath, oswMerged, params)
  fileInfo <- updateFileInfo(fileInfo, runs)
  runs <- rownames(fileInfo)
  message("Following runs will be aligned:")
  print(fileInfo[, "runName", drop=FALSE], sep = "\n")

  #### Get Precursors from the query and respectve chromatogram indices. ######
  # Get all the precursor IDs, transition IDs, Peptide IDs, Peptide Sequence Modified, Charge.
  precursors <- getPrecursors(fileInfo, oswMerged, params[["runType"]], params[["context"]],
                              params[["maxPeptideFdr"]], params[["level"]])
  if(!is.null(peps)){
    precursors <- precursors[peptide_id %in% peps, ]
    if(nrow(precursors) == 0L) stop("No peptide IDs are found in osw files.")
    setkeyv(precursors, c("peptide_id", "transition_group_id"))
  }

  #### Get OpenSWATH peak-groups and their retention times. ##########
  if(params[["transitionIntensity"]]){
    features <- getTransitions(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)
  } else{
    features <- getFeatures(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)
  }
  end_time <- Sys.time()

  #### Get the guidance tree. ####
  start_time <- Sys.time()
  distMat <- distMatrix(features, params, applyFun)
  tree <- getTree(distMat) # Check validity of tree: Names are run and master only.
  if(!is.null(newickTree)){
    tree2 <- ape::read.tree(text = newickTree)
    tree <- ape::.compressTipLabel(c(tree, tree2))[[2]] # Order tip the same way as in tree
    tree <- ape::reorder.phylo(tree, "postorder")
  }

  #### Get Peptide scores, pvalue and qvalues. ######
  masters <- tree$node.label
  peptideIDs <- unique(precursors$peptide_id)
  peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged, params[["runType"]], params[["context"]])
  peptideScores <- applyFun(peptideIDs, function(pep) {x <- peptideScores[.(pep)][,-c(1L)]
  x <- rbindlist(list(x, data.table("run" = masters, "score" = NA_real_,
                                    "pvalue" = NA_real_, "qvalue" = NA_real_)), use.names=TRUE)
  setkeyv(x, "run")
  x
  })
  names(peptideScores) <- as.character(peptideIDs)

  #### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- list2env(getMZMLpointers(fileInfo), hash = TRUE)
  message("Metadata is collected from mzML files.")

  #### Get chromatogram Indices of precursors across all runs. ############
  message("Collecting chromatogram indices for all precursors.")
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs, applyFun)
  masterChromIndex <- dummyChromIndex(precursors, masters)
  prec2chromIndex <- do.call(c, list(prec2chromIndex, masterChromIndex))

  #### Convert features into multi-peptide #####
  masterFeatures <- dummyFeatures(precursors, masters, params$transitionIntensity)
  features <- do.call(c, list(features, masterFeatures))
  start_time <- Sys.time()
  message("Building multipeptide.")
  multipeptide <- getMultipeptide(precursors, features, applyFun, NULL)
  end_time <- Sys.time()
  message("The execution time for building multipeptide:")
  print(end_time - start_time)
  message(length(multipeptide), " peptides are in the multipeptide.")

  #### Get all the child runs through hybrid alignment. ####
  adaptiveRTs <- new.env(hash = TRUE)
  refRuns <- new.env(hash = TRUE)

  # Traverse up the tree
  start_time <- Sys.time()
  traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors,
             params, adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms, applyFun)
  end_time <- Sys.time() # Report the execution time for hybrid alignment step.
  message("The execution time for creating a master run by alignment:")
  print(end_time - start_time)

  #### Map Ids from the master1 run to all parents. ####
  start_time <- Sys.time()
  traverseDown(tree, dataPath, fileInfo, multipeptide, prec2chromIndex, mzPntrs,
               precursors, adaptiveRTs, refRuns, params, applyFun)
  end_time <- Sys.time()
  message("The execution time for transfering peaks from root to runs:")
  print(end_time - start_time)

  #### Save features and add master runs to osw #####
  addMasterToOSW(dataPath, tree$node.label, oswMerged)
  save(multipeptide, fileInfo, peptideScores, refRuns, adaptiveRTs,
          file = file.path(dataPath, paste0(outFile,".temp.RData", sep = "")))

  #### Cleanup.  #######
  for(mz in names(mzPntrs)){
    if(is(mzPntrs[[mz]])[1] == "SQLiteConnection") DBI::dbDisconnect(mzPntrs[[mz]])
  }
  rm(mzPntrs, prec2chromIndex, peptideScores, features)

  #### Write tables to the disk  #######
  filename <- paste0(outFile, ".tsv", sep = "")
  finalTbl <- writeTables(fileInfo, multipeptide, precursors)
  if(params[["transitionIntensity"]]){
    finalTbl[,intensity := sapply(intensity,function(x) paste(round(x, 3), collapse=", "))]
  }
  utils::write.table(finalTbl, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Retention time alignment across runs is done.")
  message(paste0(filename, " file has been written."))

  #### End of function. #####
  alignmentStats(finalTbl, params)
  message("DONE DONE.")
}


#' Step 1 for progressive alignment
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-03-03
#' @importFrom data.table data.table setkeyv
#' @import ape
#' @inheritParams progAlignRuns
#' @seealso \code{\link{progAlignRuns}}
#' @export
progTree1 <- function(dataPath, params, outFile = "DIAlignR", oswMerged = TRUE,
                  runs = NULL, newickTree = NULL, applyFun = lapply){
  #### Check if all parameters make sense.  #########
  if(params[["chromFile"]] == "mzML"){
    if(is.null(ropenms)) stop("ropenms is required to write chrom.mzML files.")
  }
  params <- checkParams(params)

  #### Get filenames from .osw file and check consistency between osw and mzML files. #################
  fileInfo <- getRunNames(dataPath, oswMerged, params)
  fileInfo <- updateFileInfo(fileInfo, runs)
  runs <- rownames(fileInfo)
  message("Following runs will be aligned:")
  print(fileInfo[, "runName", drop=FALSE], sep = "\n")

  #### Get Precursors from the query and respectve chromatogram indices. ######
  # Get all the precursor IDs, transition IDs, Peptide IDs, Peptide Sequence Modified, Charge.
  precursors <- getPrecursors(fileInfo, oswMerged, params[["runType"]], params[["context"]],
                              params[["maxPeptideFdr"]], params[["level"]])

  #### Get OpenSWATH peak-groups and their retention times. ##########
  start_time <- Sys.time()
  if(params[["transitionIntensity"]]){
    features <- getTransitions(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)
  } else{
    features <- getFeatures(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)
  }
  message("The execution time for getting features:")
  end_time <- Sys.time()
  print(end_time - start_time)

  #### Get the guidance tree. ####
  start_time <- Sys.time()
  distMat <- distMatrix(features, params, applyFun)
  tree <- getTree(distMat) # Check validity of tree: Names are run and master only.
  if(!is.null(newickTree)){
    tree2 <- ape::read.tree(text = newickTree)
    tree <- ape::.compressTipLabel(c(tree, tree2))[[2]] # Order tip the same way as in tree
    tree <- ape::reorder.phylo(tree, "postorder")
  }

  filename <- file.path(dataPath, paste0(outFile, "_prog1.RData"))
  save(fileInfo, precursors, features, tree, file = filename, compress = FALSE)
  print(paste0("Written ", filename))

  trees <- cutTree(tree, params[["fractionNum"]])
  allNodes <- c(tree$node.label, tree$tip.label)
  mNode <- c(); rmv <- c()
  for(i in 1:params[["fractionNum"]]){
    a <- trees[[i]]
    if(is(a, "phylo")){
      mNode <- c(mNode, c(a$tip.label, a$node.label))
      tree <- ancesTree(tree, a)
    } else{
      mNode <- c(mNode, a)
    }
  }
  masters <- setdiff(allNodes, mNode)
  peptideIDs <- unique(precursors$peptide_id)
  peptideScores <- applyFun(seq_along(peptideIDs), function(i) {
    x <- data.table("run" = masters, "score" = NA_real_, "pvalue" = NA_real_, "qvalue" = NA_real_)
    setkeyv(x, "run")
    x
  })
  names(peptideScores) <- as.character(peptideIDs)
  features <- dummyFeatures(precursors, masters, params$transitionIntensity)
  prec2chromIndex <- dummyChromIndex(precursors, masters)
  multipeptide <- getMultipeptide(precursors, features, applyFun, NULL)
  outFile <- paste(outFile, 0, params[["fractionNum"]], sep = "_")
  outFile <- file.path(dataPath, paste0(outFile, ".rds"))
  tree <- ape::reorder.phylo(tree, "postorder")
  saveRDS(list(NULL, features, multipeptide, peptideScores, prec2chromIndex,
               NULL, NULL, tree), file = outFile, compress = FALSE)
  print(paste0("Written ", outFile))
  print("progTree1 is done.")
}

#' Step 2 for progressive alignment
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-03-03
#' @importFrom data.table data.table setkeyv rbindlist
#' @import DBI
#' @inheritParams progAlignRuns
#' @seealso \code{\link{progAlignRuns}}
#' @export
progSplit2 <- function(dataPath, params, outFile = "DIAlignR", oswMerged = TRUE, applyFun = lapply){
  params <- checkParams(params)
  load(file = file.path(dataPath, paste0(outFile, "_prog1.RData")))
  trees <- cutTree(tree, params[["fractionNum"]])
  tree <- trees[[params[["fraction"]]]]
  outFile <- paste(outFile, params[["fraction"]], params[["fractionNum"]], sep = "_")
  if(is(tree, "phylo")){
    fileInfo <- fileInfo[rownames(fileInfo) %in% (tree$tip.label),]
    features <- features[tree$tip.label]
    masters <- tree$node.label
  } else{
    fileInfo <- fileInfo[rownames(fileInfo) %in% tree,]
    features <- features[tree]
    masters <- NULL
  }

  #### Get Peptide scores, pvalue and qvalues. ######
  peptideIDs <- unique(precursors$peptide_id)
  peptideScores <- getPeptideScores(fileInfo, peptideIDs, oswMerged, params[["runType"]], params[["context"]])
  peptideScores <- applyFun(peptideIDs, function(pep) {x <- peptideScores[.(pep)][,-c(1L)]
  if(!is.null(masters)){
    x <- rbindlist(list(x, data.table("run" = masters, "score" = NA_real_,
                                      "pvalue" = NA_real_, "qvalue" = NA_real_)), use.names=TRUE)
  }
  setkeyv(x, "run")
  x
  })
  names(peptideScores) <- as.character(peptideIDs)

  #### Collect pointers for each mzML file. #######
  message("Collecting metadata from mzML files.")
  mzPntrs <- list2env(getMZMLpointers(fileInfo), hash = TRUE)
  message("Metadata is collected from mzML files.")

  #### Get chromatogram Indices of precursors across all runs. ############
  message("Collecting chromatogram indices for all precursors.")
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs, applyFun)
  if(!is.null(masters)){
    masterChromIndex <- dummyChromIndex(precursors, masters)
    prec2chromIndex <- do.call(c, list(prec2chromIndex, masterChromIndex))
  }

  #### Convert features into multi-peptide #####
  masterFeatures <- dummyFeatures(precursors, masters, params$transitionIntensity)
  features <- do.call(c, list(features, masterFeatures))
  start_time <- Sys.time()
  message("Building multipeptide.")
  multipeptide <- getMultipeptide(precursors, features, applyFun, NULL)
  end_time <- Sys.time()
  message("The execution time for building multipeptide:")
  print(end_time - start_time)
  message(length(multipeptide), " peptides are in the multipeptide.")

  #### Get all the child runs through hybrid alignment. ####
  adaptiveRTs <- new.env(hash = TRUE)
  refRuns <- new.env(hash = TRUE)

  #### Traverse up the tree ####
  if(!is.null(masters)){
    start_time <- Sys.time()
    traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors,
               params, adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms, applyFun)
    end_time <- Sys.time() # Report the execution time for hybrid alignment step.
    message("The execution time for creating a master run by alignment:")
    print(end_time - start_time)
    # Select the root node and save only root node data
    features <- features[masters[1]]
  }

  for(mz in names(mzPntrs)){
    if(is(mzPntrs[[mz]])[1] == "SQLiteConnection") DBI::dbDisconnect(mzPntrs[[mz]])
  }

  #### Write as rds ####
  outFile <- file.path(dataPath, paste0(outFile, ".rds"))
  saveRDS(list(fileInfo, features, multipeptide, peptideScores, prec2chromIndex,
               as.list(adaptiveRTs), as.list(refRuns), tree), file = outFile, compress = FALSE)
  print(paste0("Written ", outFile))
}

#' Step 3 for progressive alignment
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-03-03
#' @importFrom data.table setkeyv rbindlist
#' @import DBI
#' @inheritParams progAlignRuns
#' @seealso \code{\link{progAlignRuns}}
#' @export
progComb3 <- function(dataPath, params, outFile = "DIAlignR", oswMerged = TRUE, applyFun = lapply){
  params <- checkParams(params)
  #### Read all rds and RData files #####
  load(file = file.path(dataPath, paste0(outFile, "_prog1.RData")))
  peptideIDs <- unique(precursors$peptide_id)
  multipeptide <- vector(mode = "list", length = length(peptideIDs))
  peptideScores <- vector(mode = "list", length = length(peptideIDs))
  features <- prec2chromIndex <- adaptiveRTs <- refRuns <- list()
  fileInfo <- data.frame()
  for(j in 0:params[["fractionNum"]]){
    filename <- file.path(dataPath, paste0(outFile, "_", j, "_", params[["fractionNum"]], ".rds"))
    x <- readRDS(filename)
    fileInfo <- rbind(fileInfo, x[[1]])
    features <- c(features, x[[2]])
    multipeptide <- lapply(seq_along(multipeptide), function(i)
      rbindlist(list(multipeptide[[i]], x[[3]][[i]]), use.names=FALSE))
    peptideScores <- lapply(seq_along(peptideScores), function(i)
      rbindlist(list(peptideScores[[i]], x[[4]][[i]]), use.names=FALSE))
    prec2chromIndex <- c(prec2chromIndex, x[[5]])
    adaptiveRTs <- c(adaptiveRTs, x[[6]])
    refRuns <- c(refRuns, x[[7]])
    if(j == 0) tree <- x[[8]]
  }
  fileInfo <- fileInfo[tree$tip.label,]
  adaptiveRTs <- as.environment(adaptiveRTs)
  refRuns <- as.environment(refRuns)
  for (i in seq_along(peptideIDs)){
    setkeyv(multipeptide[[i]], "run")
    setkeyv(peptideScores[[i]], "run")
  }
  names(multipeptide) <- names(peptideScores) <- peptideIDs

  #### Find vertices to align ####
  message("Collecting metadata from mzML files.")
  mzPntrs <- list2env(getMZMLpointers(fileInfo), hash = TRUE)
  message("Metadata is collected from mzML files.")

  ##### Traverse to the root of all runs.  #####
  # Traverse up the tree
  start_time <- Sys.time()
  traverseUp(tree, dataPath, fileInfo, features, mzPntrs, prec2chromIndex, precursors,
             params, adaptiveRTs, refRuns, multipeptide, peptideScores, ropenms, applyFun)
  end_time <- Sys.time() # Report the execution time for hybrid alignment step.
  message("The execution time for creating a master run by alignment:")
  print(end_time - start_time)

  #### Map Ids from the master1 run to all parents. ####
  start_time <- Sys.time()
  traverseDown(tree, dataPath, fileInfo, multipeptide, prec2chromIndex, mzPntrs,
               precursors, adaptiveRTs, refRuns, params, applyFun)
  end_time <- Sys.time()
  message("The execution time for transfering peaks from root to runs:")
  print(end_time - start_time)

  for(mz in names(mzPntrs)){
    if(is(mzPntrs[[mz]])[1] == "SQLiteConnection") DBI::dbDisconnect(mzPntrs[[mz]])
  }

  #### Write data ########
  outFile <- paste(outFile, "all", params[["fractionNum"]], sep = "_")
  outFile <- file.path(dataPath, paste0(outFile, ".rds"))
  saveRDS(list(fileInfo, multipeptide, prec2chromIndex, adaptiveRTs, refRuns, features, peptideScores),
       file = outFile, compress = FALSE)
  print(paste0("Written ", outFile))
}

#' Step 4 for progressive alignment
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-03-03
#' @import DBI
#' @inheritParams progAlignRuns
#' @seealso \code{\link{progAlignRuns}}
#' @export
progSplit4 <- function(dataPath, params, outFile = "DIAlignR", oswMerged = TRUE, applyFun = lapply){
  params <- checkParams(params)
  load(file = file.path(dataPath, paste0(outFile, "_prog1.RData")))
  peptideIDs <- unique(precursors$peptide_id)
  filename <- file.path(dataPath, paste0(outFile, "_", params[["fraction"]], "_", params[["fractionNum"]], ".rds"))
  x <- readRDS(filename)
  fileInfo <- x[[1]]
  prec2chromIndex <- x[[5]]
  adaptiveRTs <-  x[[6]]
  refRuns <- x[[7]]
  tree <- x[[8]]
  message("Collecting metadata from mzML files.")
  mzPntrs <- list2env(getMZMLpointers(fileInfo), hash = TRUE)
  message("Metadata is collected from mzML files.")

  filename <- file.path(dataPath, paste0(outFile, "_", "all","_", params[["fractionNum"]], ".rds"))
  x <- readRDS(filename)
  runs <- rownames(fileInfo)
  multipeptide <- x[[2]]
  multipeptide <- lapply(seq_along(multipeptide), function(i){
    multipeptide[[i]][.(runs), ]
  })
  names(multipeptide) <- as.character(peptideIDs)

  if(is(tree, "phylo")){
    #### Map Ids from the masters to all tips. ####
    start_time <- Sys.time()
    vertices <- getNodeIDs(tree)
    ord <- rev(tree$edge[,1])
    num_merge <- length(ord)/2
    junctions <- 2*(1:num_merge)-1

    # Traverse from root to leaf node.
    for(i in junctions){
      pair <- phangorn::Descendants(tree, node = ord[i], type = "children")
      runA <- names(vertices)[vertices == pair[1]]
      runB <- names(vertices)[vertices == pair[2]]
      master <- names(vertices)[vertices == ord[i]]
      message("Mapping peaks from ", master, " to ",  runA, " and ", runB, ".")

      # Get parents to master aligned time vectors.
      filename <- file.path(dataPath, paste0(master, "_av.rds"))
      alignedVecs <- readRDS(file = filename)
      adaptiveRT <- max(adaptiveRTs[[paste(runA, runB, sep = "_")]],
                        adaptiveRTs[[paste(runB, runA, sep = "_")]])

      # Map master to runA
      refA <- refRuns[[master]][,1L][[1]]
      alignToMaster(master, runA, alignedVecs, refA, adaptiveRT, multipeptide,
                    prec2chromIndex, mzPntrs, fileInfo, precursors, params, applyFun)

      # Map master to runB
      refB <- as.integer(!(refA-1))+1L
      alignToMaster(master, runB, alignedVecs, refB, adaptiveRT, multipeptide,
                    prec2chromIndex, mzPntrs, fileInfo, precursors, params, applyFun)
    }
    end_time <- Sys.time()
    message("The execution time for transfering peaks from root to runs:")
    print(end_time - start_time)
  }

  #### Cleanup.  #######
  for(mz in names(mzPntrs)){
    if(is(mzPntrs[[mz]])[1] == "SQLiteConnection") DBI::dbDisconnect(mzPntrs[[mz]])
  }
  rm(mzPntrs, prec2chromIndex)

  #### Write tables to the disk  #######
  outFile <- paste(outFile, params[["fraction"]], params[["fractionNum"]], sep = "_")
  save(multipeptide, fileInfo, refRuns, adaptiveRTs,
       file = file.path(dataPath, paste0(outFile, ".temp.RData")))
  filename <- paste0(outFile, ".tsv", sep = "")
  finalTbl <- writeTables(fileInfo, multipeptide, precursors)
  if(params[["transitionIntensity"]]){
    finalTbl[,intensity := sapply(intensity,function(x) paste(round(x, 3), collapse=", "))]
  }
  utils::write.table(finalTbl, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Retention time alignment across runs is done.")
  message(paste0(filename, " file has been written."))

  #### End of function. #####
  alignmentStats(finalTbl, params)
  message("DONE DONE.")
}

#' Create a MST from distance matrix
#'
#' Builds a minimum spanning tree from the distance.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-05-14
#' @import ape
#' @param distMat (dist) a pairwise distance matrix.
#' @return (matrix) array of tree-edges.
#' @seealso \code{\link[ape]{mst}, \link{traverseMST}}
#' @keywords internal
#' @examples
#' m <- matrix(c(0,13,21,22, 13,0,12,13, 21,12,0,13, 22,13,13,0), byrow = TRUE,
#'             ncol = 4, dimnames = list(c("run1", "run2", "run3", "run4"),
#'                                       c("run1", "run2", "run3", "run4")))
#' distMat <- as.dist(m, diag = FALSE, upper = FALSE)
#' \dontrun{
#' x <- as.data.frame(getMST(distMat))
#' weights <- reshape2::melt(as.matrix(distMat))
#' weights$Var <- paste(weights$Var1, weights$Var2, sep = "_")
#' x$weight <- weights[match(paste(x[,1], x[,2], sep = "_"), weights$Var), "value"]
#' g1 <- igraph::graph_from_data_frame(x, directed = FALSE)
#' plot(g1, edge.label = igraph::E(g1)$weight)
#' }
getMST <- function(distMat){
  M <- ape::mst(distMat)
  M[lower.tri(M)] <- 0L
  net <- which(M == 1L, arr.ind=TRUE)
  rownames(net) <- NULL
  #runs <- attr(distMat, "Labels")
  runs <- colnames(M)
  net <- cbind("A" = runs[net[,1]], "B" = runs[net[,2]])
  net
}

#' Reorder MST
#'
#' Reorder MST for traversal from a given node.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-05-14
#' @param net (matrix) array of tree-edges.
#' @param ref0 (character) start node for traversal.
#' @return (matrix) array of tree-edges reordered for traversal.
#' @seealso \code{\link{mstAlignRuns}, \link{getMST}}
#' @keywords internal
#' @examples
#' net <- cbind("A" = c("run1","run2","run3"), "B" = c("run2","run3","run4"))
#' \dontrun{
#' traverseMST(net, "run2")
#' }
traverseMST <- function(net, ref0){
  done <- rep(FALSE, nrow(net))
  idx <- net[, 1] == ref0 | net[, 2] == ref0
  n2 <-  setdiff(c(net[idx, c(1,2)]), ref0)
  m <- cbind("ref" = ref0, "eXp" = n2)
  done[idx] <- TRUE
  while(any(done == FALSE) & length(n2) != 0L){
    ref <- n2[1]
    idx <- (net[, 1] == ref | net[, 2] == ref) & done == FALSE
    temp <-  setdiff(c(net[idx, c(1,2)]), ref)
    if(length(temp) > 0L){
      m <- rbind(m, cbind("ref" = ref, "eXp" = temp))
      done[idx] <- TRUE
    }
    n2 <- c(setdiff(n2, ref), temp)
  }
  m
}

#' Peptide quantification through MST alignment
#'
#' This function expects osw and xics directories at dataPath. It first reads osw files and fetches
#'  chromatogram indices for each analyte. To perform alignment, first a guide-tree is built using \code{\link[ape]{mst}} which
#'  can also be provided with mstNet parameter. As we traverse from the start node to the other nodes,
#'  runs are aligned pairwise.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-05-15
#' @importFrom data.table setkeyv
#' @inheritParams alignTargetedRuns
#' @param mstNet (matrix) array of tree-edges. Look up \code{\link{getMST}}.
#' @return (None)
#' @seealso \code{\link{alignTargetedRuns}, \link{progAlignRuns}, \link{getMST}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' params <- paramsDIAlignR()
#' params[["context"]] <- "experiment-wide"
#' mstAlignRuns(dataPath, params = params, outFile = "DIAlignR")
#' \dontrun{
#' y <- strsplit("run0 run1\nrun2 run2", split = '\n')[[1]]
#' y <- cbind(A = strsplit(y[1], " ")[[1]], B = strsplit(y[2], " ")[[1]])
#' plot(igraph::graph_from_edgelist(y, directed = FALSE))
#' }
#' @export
mstAlignRuns <- function(dataPath, outFile = "DIAlignR", params = paramsDIAlignR(), oswMerged = TRUE,
                         scoreFile = NULL, runs = NULL, peps = NULL, mstNet = NULL, applyFun = lapply){
  #### Check if all parameters make sense.  #########
  if(params[["chromFile"]] == "mzML"){
    if(is.null(ropenms)) stop("ropenms is required to write chrom.mzML files.")
  }
  params <- checkParams(params)

  #### Get filenames from .osw file and check consistency between osw and mzML files. #################
  fileInfo <- getRunNames(dataPath, oswMerged, params)
  fileInfo <- updateFileInfo(fileInfo, runs)
  runs <- rownames(fileInfo)
  fileInfo2 <- data.frame(fileInfo)
  if(!oswMerged) fileInfo2[["featureFile"]] <- scoreFile
  message("Following runs will be aligned:")
  print(fileInfo[, "runName", drop=FALSE], sep = "\n")

  #### Get Precursors from the query and respectve chromatogram indices. ######
  # Get all the precursor IDs, transition IDs, Peptide IDs, Peptide Sequence Modified, Charge.
  precursors <- getPrecursors(fileInfo2, oswMerged, params[["runType"]], params[["context"]],
                              params[["maxPeptideFdr"]], params[["level"]])
  if(!is.null(peps)){
    precursors <- precursors[peptide_id %in% peps, ]
    if(nrow(precursors) == 0L) stop("No peptide IDs are found in osw files.")
    setkeyv(precursors, c("peptide_id", "transition_group_id"))
  }
  if(params[["fractionNum"]] > 1L){
    idx <- getPrecursorSubset(precursors, params)
    precursors <- precursors[idx[1]:idx[2],]
    setkeyv(precursors, c("peptide_id", "transition_group_id"))
    outFile <- paste(outFile, params[["fraction"]], params[["fractionNum"]], sep = "_")
  }
  outFile <- paste0(outFile,".tsv")

  #### Get OpenSWATH peak-groups and their retention times. ##########
  start_time <- Sys.time()
  if(params[["transitionIntensity"]]){
    features <- getTransitions(fileInfo, params[["maxFdrQuery"]], params[["runType"]], applyFun)
  } else{
    features <- getFeatures(fileInfo, params[["maxFdrQuery"]], params[["maxIPFFdrQuery"]], params[["runType"]], applyFun)
  }
  end_time <- Sys.time()
  message("The execution time for fetching features:")
  print(end_time - start_time)

  #### Get the Minimum Spanning Tree. ####
  start_time <- Sys.time()
  if(is.null(mstNet)){
    distMat <- distMatrix(features, params, applyFun)
    mstNet <- getMST(distMat)
    message("Minimum spanning tree is ")
    print(paste(paste(mstNet[,1], collapse = ' '), paste(mstNet[,2], collapse = ' '), sep = '\n'))
  } else{
    y <- strsplit(mstNet, split = '\n')[[1]] # tree_count tree_rse tree_rsq tree_lin
    mstNet <- cbind(A = strsplit(y[1], " ")[[1]], B = strsplit(y[2], " ")[[1]])
  }
  nets <- lapply(runs, function(run) traverseMST(mstNet, run))
  names(nets) <- runs

  #### Get Peptide scores, pvalue and qvalues. ######
  start_time <- Sys.time()
  peptideIDs <- precursors[, logical(1), keyby = peptide_id]$peptide_id
  peptideScores <- getPeptideScores(fileInfo2, peptideIDs, oswMerged, params[["runType"]], params[["context"]])
  peptideScores <- lapply(peptideIDs, function(pep) peptideScores[.(pep)])
  names(peptideScores) <- as.character(peptideIDs)
  end_time <- Sys.time()
  message("The execution time for fetching peptide scores:")
  print(end_time - start_time)

  #### Get reference run for each precursor ########
  start_time <- Sys.time()
  message("Calculating reference run for each peptide.")
  refRuns <- getRefRun(peptideScores)
  end_time <- Sys.time()
  message("The execution time for calculating a reference run:")
  print(end_time - start_time)
  rm(peptideScores)

  #### Collect pointers for each mzML file. #######
  start_time <- Sys.time()
  message("Collecting metadata from mzML files.")
  mzPntrs <- getMZMLpointers(fileInfo)
  message("Metadata is collected from mzML files.")
  end_time <- Sys.time()
  message("The execution time for getting pointers:")
  print(end_time - start_time)

  #### Get chromatogram Indices of precursors across all runs. ############
  message("Collecting chromatogram indices for all precursors.")
  start_time <- Sys.time()
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs, applyFun)
  end_time <- Sys.time()
  message("The execution time for getting chromatogram indices:")
  print(end_time - start_time)

  #### Convert features into multi-peptide #####
  message("Building multipeptide.")
  start_time <- Sys.time()
  multipeptide <- getMultipeptide(precursors, features, params[["runType"]], applyFun, NULL)
  message(length(multipeptide), " peptides are in the multipeptide.")
  end_time <- Sys.time()
  message("The execution time for building multipeptide:")
  print(end_time - start_time)

  #### Container to save Global alignments.  #######
  message("Calculating global alignments.")
  start_time <- Sys.time()
  pairs <- lapply(1:nrow(mstNet), function(i) c(runs[mstNet[i,1]], runs[mstNet[i,2]]))

  globalFits <- lapply(1:nrow(mstNet), function(i){
    getGlobalAlignment(features, mstNet[i,1], mstNet[i,2], params[["globalAlignment"]],
                       params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
  })
  names(globalFits) <- paste(mstNet[,1], mstNet[,2], sep = "_")
  temp <- lapply(1:nrow(mstNet), function(i){
    getGlobalAlignment(features, mstNet[i,2], mstNet[i,1], params[["globalAlignment"]],
                       params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
  })
  names(temp) <- paste(mstNet[,2], mstNet[,1], sep = "_")
  globalFits <- c(globalFits, temp)
  RSE <- applyFun(globalFits, getRSE, params[["globalAlignment"]])
  globalFits <- applyFun(globalFits, extractFit, params[["globalAlignment"]])
  rm(features, temp)
  end_time <- Sys.time()
  message("The execution time for calculating global alignment:")
  print(end_time - start_time)

  #### Perform pairwise alignment ###########
  message("Performing reference-based alignment.")
  start_time <- Sys.time()
  num_of_batch <- ceiling(length(multipeptide)/params[["batchSize"]])
  invisible(
    lapply(1:num_of_batch, MSTperBatch, nets, peptideIDs, multipeptide, refRuns, precursors,
           prec2chromIndex, fileInfo, mzPntrs, params, globalFits, RSE, lapply)
  )

  #### Cleanup.  #######
  for(mz in mzPntrs){
    if(is(mz)[1] == "SQLiteConnection") DBI::dbDisconnect(mz)
    if(is(mz)[1] == "mzRpwiz") rm(mz)
  }
  rm(prec2chromIndex, globalFits, refRuns, RSE)

  end_time <- Sys.time() # Report the execution time for hybrid alignment step.
  message("The execution time for alignment:")
  print(end_time - start_time)

  #### Write tables to the disk  #######
  finalTbl <- writeTables(fileInfo, multipeptide, precursors)
  if(params[["transitionIntensity"]]){
    finalTbl[,intensity := sapply(intensity,function(x) paste(round(x, 3), collapse=", "))]
  }
  utils::write.table(finalTbl, file = outFile, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Retention time alignment across runs is done.")
  message(paste0(outFile, " file has been written."))

  #### Write alignment summary  #######
  alignmentStats(finalTbl, params)
  message("DONE DONE.")
}

#' Aligns an analyte for a batch of peptides
#'
#' For the ith analyte in the batch, this function traverse the MST from start node and align runs
#' pairwise. In the process it updates multipeptide with aligned features.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-05-15
#' @keywords internal
#' @importFrom data.table set
#' @inheritParams perBatch
#' @inherit perBatch return
#' @param nets (list) set of trees ordered for traversal obained from \code{\link{traverseMST}}.
#' @seealso \code{\link{mstAlignRuns}, \link{alignToRefMST}, \link{getAlignedTimesFast}, \link{getMultipeptide}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
MSTperBatch <- function(iBatch, nets, peptides, multipeptide, refRuns, precursors, prec2chromIndex,
                     fileInfo, mzPntrs, params, globalFits, RSE, applyFun = lapply){
  # if(params[["chromFile"]] =="mzML") fetchXIC = extractXIC_group
  fetchXICs = extractXIC_group2
  message("Processing Batch ", iBatch)
  batchSize <- params[["batchSize"]]
  strt <- ((iBatch-1)*batchSize+1)
  stp <- min((iBatch*batchSize), length(peptides))
  runs <- rownames(fileInfo)

  ##### Get XICs into memory for the batch across all runs #####
  pIdx <- lapply(peptides[strt:stp], function(pep) which(precursors$peptide_id == pep))
  analytesA <- lapply(pIdx, function(i) .subset2(precursors, "transition_group_id")[i])
  chromIndices <- lapply(runs, function(run) lapply(pIdx, function(i) .subset2(prec2chromIndex[[run]], "chromatogramIndex")[i]))
  cons <- lapply(seq_along(runs), function(i) createTemp(mzPntrs[[runs[i]]], unlist(chromIndices[[i]])))
  names(cons) <- names(chromIndices) <- runs

  ##### Get aligned multipeptide for the batch #####
  invisible(applyFun(strt:stp, function(rownum){
    peptide <- peptides[rownum]
    DT <- multipeptide[[rownum]]
    ref <- refRuns[rownum, "run"][[1]]
    idx <- (rownum - (iBatch-1)*batchSize)
    analytes <- analytesA[[idx]]
    net <- nets[[ref]]

    XICs <- lapply(seq_along(runs), function(i){
      cI <- chromIndices[[i]][[idx]]
      if(any(is.na(unlist(cI))) | is.null(unlist(cI))) return(NULL)
      temp <- lapply(cI, function(i1) fetchXICs(cons[[i]], i1))
      names(temp) <- as.character(analytes)
      temp
    })
    names(XICs) <- runs

    XICs.ref <- XICs[[ref]]
    if(is.null(XICs.ref)){
      message("Chromatogram indices for peptide ", peptide, " are missing in ", fileInfo[ref, "runName"])
      message("Skipping peptide ", peptide, " across all runs.")
      return(invisible(NULL))
      # Select refrence from net
    }
    ##### Set alignment rank for all precrusors of the peptide in the reference run #####
    refIdx <- which(DT[["run"]] == ref & DT[["peak_group_rank"]] == 1L)
    refIdx <- refIdx[which.min(DT$m_score[refIdx])]
    if(length(refIdx)==0) {
      message("Features for peptide ", peptide, " is missing in ", fileInfo[ref, "runName"])
      message("Skipping peptide ", peptide, " across all runs.")
      return(invisible(NULL))
    }
    set(DT, i = refIdx, 10L, 1L)
    setOtherPrecursors(DT, refIdx, XICs.ref, analytes, params)

    ##### Align all runs to reference run and set their alignment rank #####
    invisible(
      lapply(1:nrow(net), alignToRefMST, net, fileInfo, XICs, params, analytes,
             DT, globalFits, RSE)
    )

    ##### Return the dataframe with alignment rank set to TRUE #####
    updateOnalignTargetedRuns(rownum)
  })
  )
  for(con in cons) DBI::dbDisconnect(con)
  invisible(NULL)
}

#' Aligns an analyte for an edge of MST
#'
#' df contains unaligned features for an analyte across multiple runs. This function aligns eXp run to
#' ref run and updates corresponding features.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2021) + GPL-3
#' Date: 2021-05-15
#' @keywords internal
#' @inheritParams alignToRef
#' @inherit alignToRef return
#' @param iNet (integer) the index of edge to be aligned in the net.
#' @param net (matrix) each row represents an edge of MST.
#' @param analytes (string) precursor IDs of the requested peptide.
#' @param df (dataframe) a collection of features related to analytes.
#' @seealso \code{\link{mstAlignRuns}, \link{MSTperBatch}, \link{setAlignmentRank}, \link{getMultipeptide}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
alignToRefMST <- function(iNet, net, fileInfo, XICs, params, analytes,
                          df, globalFits, RSE){
  ref <- net[[iNet, 1]]; eXp <- net[[iNet, 2]]
  # Get XIC_group from experiment run.
  XICs.eXp <- XICs[[eXp]]
  eXpIdx <- which(df[["run"]] == eXp)
  ##### Check if any feature is below unaligned FDR. If present alignment_rank = 1. #####
  if(any(.subset2(df, "m_score")[eXpIdx] <=  params[["unalignedFDR"]], na.rm = TRUE)){
    tempi <- eXpIdx[which.min(df$m_score[eXpIdx])]
    set(df, tempi, 10L, 1L)
    if(is.null(XICs.eXp)) return(invisible(NULL))
    setOtherPrecursors(df, tempi, XICs.eXp, analytes, params)
    return(invisible(NULL))
  }

  # No high quality feature, hence, alignment is needed.
  # check is alignment rank is set in reference.
  refIdx <- which(df[["run"]] == ref)
  refIdx <- refIdx[which(.subset2(df, 10L)[refIdx] == 1L)]
  if(length(refIdx) == 0L) return(invisible(NULL))
  ss <- .subset2(df, "m_score")[refIdx]
  refIdx <- ifelse(all(is.na(ss)), refIdx[1], refIdx[which.min(ss)])

  ### if XICs are missing, go to next run. ####
  XICs.ref <- XICs[[ref]]
  if(is.null(XICs.ref) || is.null(XICs.eXp)){
    message("Chromatogram indices for precursor ", analytes, " are missing in either ",
            fileInfo[ref, "runName"], " or", fileInfo[eXp, "runName"])
    message("Skipping precursor ", analytes, " in the alignment of this pair.")
    return(invisible(NULL))
  }

  # Select 1) all precursors OR 2) high quality precursor
  if(FALSE){
    # Turned off as precursor XICs have different time ranges.
    XICs.ref.pep <- unlist(XICs.ref, recursive = FALSE, use.names = FALSE)
    XICs.eXp.pep <- unlist(XICs.eXp, recursive = FALSE, use.names = FALSE)
  } else {
    analyte_chr <- as.character(.subset2(df, 1L)[[refIdx]])
    XICs.ref.pep <- XICs.ref[[analyte_chr]]
    XICs.eXp.pep <- XICs.eXp[[analyte_chr]]
  }

  ##### Get the aligned Indices #####
  pair <- paste(ref, eXp, sep = "_")
  globalFit <- globalFits[[pair]]
  adaptiveRT <- params[["RSEdistFactor"]]*RSE[[pair]]

  if(missingInXIC(XICs.eXp.pep) || missingInXIC(XICs.ref.pep)){
    message("Missing values in the chromatogram of ", paste0(analytes, sep = " "), "precursors in run either ",
            fileInfo[ref, "runName"], " or", fileInfo[eXp, "runName"])
    return(invisible(NULL)) # Missing values in chromatogram
  }

  tAligned <- tryCatch(expr = getAlignedTimesFast(XICs.ref.pep, XICs.eXp.pep, globalFit, adaptiveRT,
                                                  params),
                       error = function(e){
                         message("\nError in the alignment of ", paste0(analytes, sep = " "), "precursors in runs ",
                                 fileInfo[ref, "runName"], " and ", fileInfo[eXp, "runName"])
                         warning(e)
                         return(NULL)
                       })
  if(is.null(tAligned)) return(invisible(NULL))
  tryCatch(expr = setAlignmentRank(df, refIdx, eXp, tAligned, XICs.eXp, params, adaptiveRT),
           error = function(e){
             message("\nError in setting alignment rank of ", paste0(analytes, sep = " "), "precursors in runs ",
                     fileInfo[eXp, "runName"], " and ", fileInfo[eXp, "runName"])
             warning(e)
             return(invisible(NULL))
           })

  tempi <- eXpIdx[which(df$alignment_rank[eXpIdx] == 1L)]
  if(length(tempi) == 0L) return(invisible(NULL))
  setOtherPrecursors(df, tempi, XICs.eXp, analytes, params)
  invisible(NULL)
}

#!/bin/usr/Rscript

#' Command line interface to alignTargetedRuns
#'
#' This function is to be used via the command line, to run alignTargetedRuns
#' @author Justin Sing, \email{justinc.sing@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-0386-0092
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2021-08-27
#' @inheritParams alignTargetedRuns
#' @seealso \code{\link{alignTargetedRuns}
#' @examples
#' Rscript alignTargetedRuns_cli.R --dataPath=/data/osw/ --params=context:experiment-wide,maxFdrQuery:0.01
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
alignTargetedRuns_cli <- function(){
  # Get Arguments passed to command line
  args <- commandArgs(trailingOnly = TRUE)
  # Extract Passed Command Line Arguments
  args <- parseArgs(args)
  # Register BioCParallel workers if being used
  if (("regBioCP" %in% names(args))){
    message("INFO: Registering BiocParallel...")
    eval(parse(text = args[["regBioCP"]]))
  }
  # Main call to DIAlignR's alingment of Targeted Runs
  # alignTargetedRuns(dataPath=args[["dataPath"]], outFile=args[["outFile"]], params=args[["params"]], oswMerged=args[["oswMerged"]],
  #                   runs=args[["runs"]], peps=args[["peps"]], refRun=args[["refRun"]], applyFun=args[["applyFun"]])
}

#' Parse command line arguments and create named argument list
#'
#' @author Justin Sing, \email{justinc.sing@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-0386-0092
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2021-08-27
#' @importFrom BiocParallel bplapply register MulticoreParam
#' @inheritParams paramsDIAlignR
#' @seealso \code{\link{paramsDIAlignR}
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
parseArgs <- function(args){
  # "--dataPath=/data/stuff/ --oswMerged=TRUE --run=run0,run1,run2 --pep=PEP,TIDE --refRun=run4 --applyFun=lapply"
  # args <- c("--dataPath=./", "--outFile=dialignr", "--oswMerged=TRUE", "--runs=run0,run1,run2", "--peps=0,1", "--refRun=run4",
  #           "--applyFun=BiocParallel::bplapply", "--regBioCP=BiocParallel::register(BiocParallel::MulticoreParam(workers=4,progressbar=TRUE))",
  #           "--params=context:experiment-wide,maxFdrQuery:0.01,fitEMG:TRUE")
  args_list <- strsplit(args, "\\s")
  args_list <- gsub("--", "", args_list)

  # Display help message
  if ( "help" %in% args_list ){
    display_help()
    stop(call. = FALSE)
  }

  # Check for unknown args, and ignore them if present
  unknown_args <- args_list[!grepl("dataPath|outFile|oswMerged|runs|peps|refRun|applyFun|regBioCP|params", args_list)]
  if (length(unknown_args)>0){
    warning(sprintf("WARN: There was an unknown argument passed (--%s), going to ignore it!\n", unknown_args))
    # Remove unknown args
    args_list <- args_list[ !(args_list %in% unknown_args) ]
  }

  args_keys <- unlist(lapply(args_list, function(x) { regmatches(x, regexpr("=", x), invert = TRUE)[[1]][[1]] }))
  args_values <- lapply(args_list, function(x) { regmatches(x, regexpr("=", x), invert = TRUE)[[1]][[2]] })
  names(args_values) <- args_keys

  # Ensure dataPath is given
  stopifnot("--dataPath was not provided!" = any(args_keys %in% c("dataPath")))

  # Check Arguments
  ## Check to make sure dataPath is valid
  stopifnot("Provided dataPath does not exist!" = file.exists( path.expand(normalizePath(args_values[["dataPath"]])) ) )
  args_values[["dataPath"]] <- path.expand(normalizePath(args_values[["dataPath"]]))
  ## Check to see if outFile is in list
  if ( !("outFile" %in% args_keys) ){
    message("INFO: Setting outFile name to default dialignr")
    args_values[["outFile"]] <- "dialignr"
  }
  ## Check to see if oswMerged is in list
  if ( "oswMerged" %in% args_keys ){
    args_values[["oswMerged"]] <- as.logical(args_values[["oswMerged"]])
    ## Check to make sure oswMerged is valid
    stopifnot("oswMerged can only be TRUE or FALSE!" = isTRUE(args_values[["oswMerged"]]) | isFALSE(args_values[["oswMerged"]]) )
  } else {
    message("INFO: Setting oswMerged to default TRUE")
    args_values[["oswMerged"]] <- TRUE
  }
  ## Check to see if run is in list
  if ( "runs" %in% args_keys ){
    ## Convert run to character vector
    args_values[["runs"]] <- unlist(strsplit(args_values[["runs"]], ","))
  } else {
    message("INFO: Setting runs to default NULL")
    args_values[["runs"]] <- NULL
  }
  ## Check to see if refRun is in list
  if ( !("refRun" %in% args_keys) ){
    message("INFO: Setting refRun to default NULL")
    args_values[["refRun"]] <- NULL
  }
  ## Check to see if peps is in list
  if ( "peps" %in% args_keys ){
    ## Convert pep to character vector
    args_values[["peps"]] <- unlist(lapply(unlist(strsplit(args_values[["peps"]], ",")), check_and_convert_class))
  } else {
    message("INFO: Setting peps to default NULL")
    args_values[["peps"]] <- NULL
  }
  ## Check to see if applyFun is in list
  if ( "applyFun" %in% args_keys ){
    if ( ("BiocParallel::bplapply" %in% args_values[["applyFun"]]) | ("lapply" %in% args_values[["applyFun"]]) ){
      if ( ("BiocParallel::bplapply" %in% args_values[["applyFun"]]) & !("regBioCP" %in% args_keys)){
        warning(sprintf("WARN: --regBioCP was not set, will use register cores using this: BiocParallel::register(BiocParallel::MulticoreParam(workers=4,progressbar=TRUE))"))
        args_values[["regBioCP"]] <- "BiocParallel::register(BiocParallel::MulticoreParam(workers=4,progressbar=TRUE))"
      }
      args_values[["applyFun"]] <- eval(parse(text = args_values[["applyFun"]]))
    } else {
      warning(sprintf("WARN: Supplied --applyFun=%s, has to be lappy or BiocParallel::bplapply, setting to default lapply", args_values[["applyFun"]]))
      args_values[["applyFun"]] <- lapply
    }

  } else {
    message("INFO: Setting applyFun to default lapply")
    args_values[["applyFun"]] <- lapply
  }
  ## Check to see if params is in list
  if ( "params" %in% args_keys ){
    # Get list of default params to update
    params <- DIAlignR::paramsDIAlignR()
    # Get params to change
    new_params <- unlist(strsplit(args_values[["params"]], ","))
    new_params_keys <- unlist(lapply(new_params, function(x) { strsplit(x, ":")[[1]][[1]] }))
    new_params_values <- lapply(new_params, function(x) { check_and_convert_class(strsplit(x, ":")[[1]][[2]]) })
    names(new_params_values) <- new_params_keys
    # Update params with new values
    for (new_params_key in new_params_keys){
      message(sprintf("INFO: Updating params[['%s']] from %s -> %s", new_params_key, params[[new_params_key]], new_params_values[[new_params_key]]))
      params[[new_params_key]] <- new_params_values[[new_params_key]]
    }
    args_values[["params"]] <- params
  } else {
    message("INFO: Setting params to default DIAlignR::paramsDIAlignR()")
    args_values[["params"]] <- DIAlignR::paramsDIAlignR()
  }
  return(args_values)
}

#' Convert a character vector to logical, numeric, or stay as character
#'
#' @author Justin Sing, \email{justinc.sing@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-0386-0092
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2021-08-27
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
check_and_convert_class <- function(x){
  suppressWarnings({
    if (is.logical(as.logical(x)) & !is.na(as.logical(x))) {
      return(as.logical(x))
    } else if (is.numeric(as.numeric(x)) & !is.na(as.numeric(x))) {
      return(as.numeric(x))
    } else {
      return(x)
    }
  })
}

#' Display help message for running CLI
#'
#' @author Justin Sing, \email{justinc.sing@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-0386-0092
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2021-08-27
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
display_help <- function(){
  message(
    "Name: 
        Run DIAlignR's alignTargetedRuns via the Command Line

    Usage: 
        Rscript alignTargetedRuns_cli.R --dataPath=/data/ [args] | --help

        Example: Rscript alignTargetedRuns_cli.R --dataPath=/data/osw/ --params=context:experiment-wide,maxFdrQuery:0.01

        Example2: Rscript alignTargetedRuns_cli.R --dataPath=/data/osw/ --oswMerged=FALSE --params=context:experiment-wide,maxFdrQuery:0.01 --runs=run0,run1,run2 --peps=0,1 --applyFun=BiocParallel::bplapply --regBioCP=BiocParallel::register(BiocParallel::MulticoreParam(workers=4,progressbar=TRUE))

    Options: 
        --dataPath: path to xics and osw directory.
        --outFile: name of the output file.
        --oswMerged: TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
        --params: Parameters for the alignment functions generated from DIAlignR::paramsDIAlignR(). Separate keys and values using a ':', and separate parameters using ','. Example: --params=context:experiment-wide,maxFdrQuery:0.01,fitEMG:TRUE
        --runs: names of xics file without extension. Separate runs using ','. Example: --runs=run0,run1,run2
        --refRun: reference for alignment. If no run is provided, m-score is used to select reference run.
        --peps: ids of peptides to be aligned. If NULL, align all peptides. Separate peptide ids using ','. Example--peps=1,2,3
        --appyFun: value must be either lapply or BiocParallel::bplapply.
        --regBioCP: If using BiocParallel::bplapply, register cores to use. Example: --regBioCP=BiocParallel::register(BiocParallel::MulticoreParam(workers=4,progressbar=TRUE)) . Make sure there are no spaces in this command
        --help: Display this help message
    "
  )
}

#' Main CLI Method Call
#'
#' @author Justin Sing, \email{justinc.sing@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-0386-0092
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2021-08-27
#' @inheritParams alignTargetedRuns
#' @seealso \code{\link{alignTargetedRuns}
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
#' @export
alignTargetedRuns_cli()

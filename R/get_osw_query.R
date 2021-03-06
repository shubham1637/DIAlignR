#' Generate SQL query to fetch information from osw files.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param maxFdrQuery (numeric) value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param analytes (vector of strings) transition_group_ids for which features are to be extracted. analyteInGroupLabel must be set according the pattern used here.
#' @param filename (string) as mentioned in RUN table of osw files.
#' @param runType (char) This must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @param identifying (logical) TRUE for the extraction of identifying transtions. (Default: FALSE)
#' @return SQL query to be searched.
#' @keywords internal
getQuery <- function(maxFdrQuery, oswMerged = TRUE, analytes = NULL,
                     filename = NULL, runType = "DIA_Proteomics", analyteInGroupLabel = FALSE, identifying=FALSE){
  if(is.null(analytes)){
    selectAnalytes <- ""
  } else{
    selectAnalytes <- paste0(" AND transition_group_id IN ('", paste(analytes, collapse="','"),"')")
  }

  if(oswMerged){
    matchFilename <- paste0(" AND RUN.FILENAME ='", filename,"'")
  } else{
    matchFilename <- ""
  }

  if(analyteInGroupLabel == TRUE){
    transition_group_id <- " PRECURSOR.GROUP_LABEL AS transition_group_id"
  } else {
    transition_group_id <- " PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.CHARGE AS transition_group_id"
  }

  if(runType == "DIA_Metabolomics"){
    query <- paste0("SELECT RUN.ID AS id_run,
    COMPOUND.ID AS id_compound,
    COMPOUND.COMPOUND_NAME || '_' || COMPOUND.ADDUCTS AS transition_group_id,
    TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
    RUN.ID AS run_id,
    RUN.FILENAME AS filename,
    FEATURE.EXP_RT AS RT,
    FEATURE.DELTA_RT AS delta_rt,
    PRECURSOR.LIBRARY_RT AS assay_RT,
    FEATURE.ID AS id,
    COMPOUND.SUM_FORMULA AS sum_formula,
    COMPOUND.COMPOUND_NAME AS compound_name,
    COMPOUND.ADDUCTS AS Adducts,
    PRECURSOR.CHARGE AS Charge,
    PRECURSOR.PRECURSOR_MZ AS mz,
    FEATURE_MS2.AREA_INTENSITY AS Intensity,
    FEATURE.LEFT_WIDTH AS leftWidth,
    FEATURE.RIGHT_WIDTH AS rightWidth,
    SCORE_MS2.RANK AS peak_group_rank,
    SCORE_MS2.QVALUE AS m_score
    FROM PRECURSOR
    INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID
    INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID
    INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
    INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
    LEFT JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
    LEFT JOIN FEATURE_MS1 ON FEATURE_MS1.FEATURE_ID = FEATURE.ID
    LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
    LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
    WHERE COMPOUND.DECOY = 0 AND SCORE_MS2.QVALUE <  ", maxFdrQuery, selectAnalytes, matchFilename, "
    ORDER BY transition_group_id,
    peak_group_rank;")
  } else if (runType == "MRM_Proteomics"){
    query <- paste0("SELECT PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.CHARGE AS transition_group_id,
  RUN.FILENAME AS filename,
  FEATURE.EXP_RT AS RT,
  FEATURE.DELTA_RT AS delta_rt,
  PRECURSOR.LIBRARY_RT AS assay_RT,
  FEATURE_MS2.AREA_INTENSITY AS Intensity,
  FEATURE.LEFT_WIDTH AS leftWidth,
  FEATURE.RIGHT_WIDTH AS rightWidth,
  TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id
  FROM PRECURSOR
  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AND PRECURSOR.DECOY=0
  INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
  LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
  ORDER BY transition_group_id;")
  } else if (runType=="DIA_IPF"){
    join_score_ms2 <- "INNER JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID"
    join_score_ipf <- sprintf("INNER JOIN SCORE_IPF ON SCORE_IPF.FEATURE_ID = FEATURE.ID")
    join_peptide <- sprintf("INNER JOIN PEPTIDE ON ( PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID OR PEPTIDE.ID = SCORE_IPF.PEPTIDE_ID )")
    miscellaneous_score_ipf_control <- "AND SCORE_IPF.QVALUE IS NOT NULL AND PEPTIDE.MODIFIED_SEQUENCE NOT LIKE '%UniMod%'"
    if ( identifying ){
      ## Filter Identifying transitions for PEP level threshold, and keep detecting NULL transitions
      identifying_transition_filter_query <- sprintf("AND (SCORE_TRANSITION.PEP < %s OR (TRANSITION.DETECTING AND SCORE_TRANSITION.PEP IS NULL))", identifying.transitionPEPfilter)
    } else {
      identifying_transition_filter_query <- ''
    }
    query <- paste0("SELECT", transition_group_id, ",
        RUN.FILENAME AS filename,
        FEATURE.ID AS feature_id,
        FEATURE.EXP_RT AS RT,
        FEATURE.DELTA_RT AS delta_rt,
        FEATURE.EXP_RT - FEATURE.DELTA_RT AS assay_rt,
        FEATURE_MS2.AREA_INTENSITY AS Intensity,
        FEATURE.LEFT_WIDTH AS leftWidth,
        FEATURE.RIGHT_WIDTH AS rightWidth,
        SCORE_MS2.RANK AS peak_group_rank,
        SCORE_MS2.SCORE AS d_score,
        SCORE_MS2.PEP AS ms2_pep,
        SCORE_MS2.QVALUE as ms2_m_score,
        SCORE_IPF.PEP AS ipf_pep,
        SCORE_IPF.QVALUE AS m_score,
        TRANSITION.ID AS transition_id,
        TRANSITION.PRODUCT_MZ AS product_mz,
        ---SCORE_TRANSITION.FEATURE_ID AS score_transition_feature_id,
        ---SCORE_TRANSITION.TRANSITION_ID AS score_transition_id,
        ---SCORE_TRANSITION.PEP AS transition_pep,
        TRANSITION.DETECTING AS detecting_transitions,
        TRANSITION.IDENTIFYING AS identifying_transitions
  FROM FEATURE
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  INNER JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID\n",
  join_score_ms2,
  "\n",
  join_score_ipf,
  "\nINNER JOIN PRECURSOR ON PRECURSOR.ID = FEATURE.PRECURSOR_ID
  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = FEATURE.PRECURSOR_ID\n",
  join_peptide,
  "\nINNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN TRANSITION ON TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID = TRANSITION.ID
  WHERE FEATURE.ID IS NOT NULL -- Default WHERE Being clause
  AND SCORE_IPF.QVALUE < ", maxFdrQuery,
  "\n",
  miscellaneous_score_ipf_control,
  "\n",
  identifying_transition_filter_query,
  "\n",
  selectAnalytes,
  "\n",
  matchFilename,
  "\nAND (
        TRANSITION.DETECTING=TRUE
        OR TRANSITION.IDENTIFYING=", identifying,
            ")
  ORDER BY transition_group_id,
  peak_group_rank
  ")
  }else{
    query <- paste0("SELECT", transition_group_id,",
  RUN.FILENAME AS filename,
  FEATURE.EXP_RT AS RT,
  FEATURE.DELTA_RT AS delta_rt,
  PRECURSOR.LIBRARY_RT AS assay_RT,
  FEATURE_MS2.AREA_INTENSITY AS Intensity,
  FEATURE.LEFT_WIDTH AS leftWidth,
  FEATURE.RIGHT_WIDTH AS rightWidth,
  SCORE_MS2.RANK AS peak_group_rank,
  SCORE_MS2.QVALUE AS m_score,
  TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id
  FROM PRECURSOR
  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AND PRECURSOR.DECOY=0
  INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
  LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
  LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
  WHERE SCORE_MS2.QVALUE < ", maxFdrQuery, selectAnalytes, matchFilename, "
  ORDER BY transition_group_id,
  peak_group_rank;")
  }
  return(query)
}


#' Generate SQL query to fetch limited information from osw files.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-14
#' @param maxFdrQuery (numeric) value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param filename (string) as mentioned in RUN table of osw files..
#' @param runType (char) This must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @return SQL query to be searched.
#' @seealso \code{\link{getOswAnalytes}}
#' @keywords internal
getAnalytesQuery <- function(maxFdrQuery, oswMerged = TRUE, filename = NULL,
                             runType = "DIA_Proteomics", analyteInGroupLabel = FALSE){
  if(oswMerged){
    matchFilename <- paste0(" AND RUN.FILENAME ='", filename,"'")
  } else{
    matchFilename <- ""
  }

  if(analyteInGroupLabel == TRUE){
    transition_group_id <- " PRECURSOR.GROUP_LABEL AS transition_group_id"
  } else {
    transition_group_id <- " PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.CHARGE AS transition_group_id"
  }

  if(runType == "DIA_Metabolomics"){
    query <- paste0("SELECT COMPOUND.ID AS compound_id,
    COMPOUND.COMPOUND_NAME || '_' || COMPOUND.ADDUCTS AS transition_group_id,
    RUN.FILENAME AS filename,
    SCORE_MS2.RANK AS peak_group_rank,
    SCORE_MS2.QVALUE AS m_score
    FROM PRECURSOR
    INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID
    INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID
    INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
    INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
    LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
    WHERE COMPOUND.DECOY = 0 AND SCORE_MS2.QVALUE <  ", maxFdrQuery, matchFilename, "
    ORDER BY transition_group_id,
    peak_group_rank;")
  } else if (runType == "MRM_Proteomics"){
    query <- paste0("SELECT PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.CHARGE AS transition_group_id,
  RUN.FILENAME AS filename,
  FEATURE.EXP_RT AS RT,
  FROM PRECURSOR
  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AND PRECURSOR.DECOY=0
  INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  ORDER BY transition_group_id;")
  } else{
    query <- paste0("SELECT", transition_group_id,",
  RUN.FILENAME AS filename,
  SCORE_MS2.RANK AS peak_group_rank,
  SCORE_MS2.QVALUE AS m_score,
  TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id
  FROM PRECURSOR
  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AND PRECURSOR.DECOY=0
  INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
  LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
  WHERE SCORE_MS2.QVALUE < ", maxFdrQuery, matchFilename, "
  ORDER BY transition_group_id,
  peak_group_rank;")
  }
  return(query)
}

#  https://stackoverflow.com/questions/10622260/how-do-you-query-an-int-column-for-any-value
#' Get precursor Info
#'
#' For each precursor in the table respective transition ids are fetched.
#' Order of transition is kept same as the order of their intensities in \code{\link{getTransitionsQuery}}.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2020-04-04
#' @inheritParams fetchPrecursorsInfo
#' @return SQL query to be searched.
#' @seealso \code{\link{fetchPrecursorsInfo}}
#' @keywords internal
getPrecursorsQuery <- function(runType = "DIA_Proteomics", level = "Peptide")
{
  if (runType == "DIA_Proteomics"){
    if(level == "Protein"){
      query <- "SELECT DISTINCT PRECURSOR.ID AS transition_group_id,
      TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
      PEPTIDE.ID AS peptide_id,
      PEPTIDE.MODIFIED_SEQUENCE AS sequence,
      PRECURSOR.CHARGE AS charge,
      PRECURSOR.GROUP_LABEL AS group_label
      FROM PEPTIDE_PROTEIN_MAPPING
      INNER JOIN (
            SELECT PROTEIN_ID
            FROM SCORE_PROTEIN
            WHERE SCORE_PROTEIN.CONTEXT = $CONTEXT AND SCORE_PROTEIN.QVALUE < $FDR
            ) AS SCORE_PROTEIN ON SCORE_PROTEIN.PROTEIN_ID = PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID
       INNER JOIN(
            SELECT PEPTIDE_ID
            FROM SCORE_PEPTIDE
            WHERE SCORE_PEPTIDE.CONTEXT = $CONTEXT AND SCORE_PEPTIDE.QVALUE < $FDR
            ) AS SCORE_PEPTIDE ON SCORE_PEPTIDE.PEPTIDE_ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID
      INNER JOIN(
      	  SELECT ID, MODIFIED_SEQUENCE
      	  FROM PEPTIDE
      	  WHERE PEPTIDE.DECOY = 0
      	  ) AS PEPTIDE ON PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID = PEPTIDE.ID
      LEFT JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
      INNER JOIN (
      		SELECT ID, CHARGE, GROUP_LABEL
      		FROM PRECURSOR
      		WHERE PRECURSOR.DECOY = 0
      	  ) AS PRECURSOR ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      LEFT JOIN TRANSITION_PRECURSOR_MAPPING  ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      ORDER BY peptide_id, transition_group_id, transition_id;"
    } else{
      query <- "SELECT DISTINCT PRECURSOR.ID AS transition_group_id,
      TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
      PEPTIDE.ID AS peptide_id,
      PEPTIDE.MODIFIED_SEQUENCE AS sequence,
      PRECURSOR.CHARGE AS charge,
      PRECURSOR.GROUP_LABEL AS group_label
      FROM PRECURSOR
      INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
      INNER JOIN (
      SELECT PEPTIDE_ID
      FROM SCORE_PEPTIDE
      WHERE SCORE_PEPTIDE.CONTEXT = $CONTEXT AND SCORE_PEPTIDE.QVALUE < $FDR
      ) AS SCORE_PEPTIDE ON SCORE_PEPTIDE.PEPTIDE_ID = PEPTIDE.ID
      WHERE PRECURSOR.DECOY = 0
      ORDER BY peptide_id, transition_group_id, transition_id;"
    }
  } else if (runType == "DIA_IPF"){
    ## TODO: Not as clear on PROTEIN level, for now we only do peptide level
      query <- "SELECT DISTINCT PRECURSOR.ID AS transition_group_id,
      TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
      PEPTIDE.ID AS peptide_id,
      PEPTIDE.MODIFIED_SEQUENCE AS sequence,
      PRECURSOR.CHARGE AS charge,
      PRECURSOR.GROUP_LABEL AS group_label
      --TRANSITION.DETECTING AS detecting,
	    --TRANSITION.IDENTIFYING AS identifying
      FROM PRECURSOR
      INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN (
      SELECT *
      FROM TRANSITION
      WHERE (
      TRANSITION.DETECTING=TRUE
      OR TRANSITION.IDENTIFYING=$USE_IDENTIFYING --- #identifying
          )
      ) AS TRANSITION ON TRANSITION.ID = TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID
      INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
      INNER JOIN (
      SELECT PEPTIDE_ID
      FROM SCORE_PEPTIDE
      WHERE SCORE_PEPTIDE.CONTEXT = $CONTEXT AND SCORE_PEPTIDE.QVALUE < $FDR
      ) AS SCORE_PEPTIDE ON SCORE_PEPTIDE.PEPTIDE_ID = PEPTIDE.ID
      WHERE PRECURSOR.DECOY = 0
      ORDER BY peptide_id, transition_group_id, transition_id;"
  } else if (runType == "DIA_Metabolomics") {
      query <- "SELECT PRECURSOR.ID AS transition_group_id,
      TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
      PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID AS peptide_id,
      COMPOUND.SUM_FORMULA AS sum_formula,
      COMPOUND.COMPOUND_NAME AS compound_name,
      COMPOUND.ADDUCTS AS adducts,
      PRECURSOR.CHARGE AS charge,
      PRECURSOR.GROUP_LABEL AS group_label
      FROM PRECURSOR
      INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID = PRECURSOR.ID
      INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID
      WHERE PRECURSOR.DECOY = 0
      ORDER BY peptide_id, transition_group_id, transition_id;"
  }
  query
}

#' Get features from a SQLite file
#'
#' Query is generated to identify features below a FDR cut-off from a run.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2020-04-07
#' @param runType (char) This must be one of the strings "DIA_Proteomics", "DIA_IPF", "DIA_Metabolomics".
#' @return SQL query to be searched.
#' @seealso \code{\link{fetchFeaturesFromRun}}
#' @keywords internal
getFeaturesQuery <- function(runType = "DIA_Proteomics") #is the same call for both
{
  if (runType == "DIA_Proteomics"){
    query <- "SELECT PRECURSOR.ID AS transition_group_id,
    FEATURE.ID AS feature_id,
    FEATURE.EXP_RT AS RT,
    FEATURE_MS2.AREA_INTENSITY AS intensity,
    FEATURE.LEFT_WIDTH AS leftWidth,
    FEATURE.RIGHT_WIDTH AS rightWidth,
    SCORE_MS2.RANK AS peak_group_rank,
    SCORE_MS2.QVALUE AS m_score
    FROM PRECURSOR
    INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
    INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
    LEFT JOIN (
        SELECT FEATURE_ID, AREA_INTENSITY
        FROM FEATURE_MS2
    ) AS FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
    INNER JOIN (
        SELECT FEATURE_ID, RANK, QVALUE
        FROM SCORE_MS2
        WHERE SCORE_MS2.QVALUE < $FDR
        ) AS SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
    WHERE PRECURSOR.DECOY = 0 AND RUN.ID = $runID
    ORDER BY transition_group_id, peak_group_rank;"
  } else if (runType == "DIA_IPF"){
    query <- "SELECT PRECURSOR.ID AS transition_group_id,
    FEATURE.ID AS feature_id,
    FEATURE.EXP_RT AS RT,
    FEATURE_MS2.AREA_INTENSITY AS intensity,
    FEATURE.LEFT_WIDTH AS leftWidth,
    FEATURE.RIGHT_WIDTH AS rightWidth,
    SCORE_MS2.RANK AS peak_group_rank,
    SCORE_MS2.QVALUE AS ms2_m_score,
    SCORE_IPF.QVALUE AS m_score
    FROM PRECURSOR
    INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
    INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
    INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
    LEFT JOIN (
        SELECT FEATURE_ID, AREA_INTENSITY
        FROM FEATURE_MS2
    ) AS FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
    INNER JOIN (
        SELECT FEATURE_ID, RANK, QVALUE
        FROM SCORE_MS2
        WHERE SCORE_MS2.QVALUE < $FDR
        ) AS SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
    INNER JOIN (
        SELECT FEATURE_ID, QVALUE, PEPTIDE_ID
        FROM SCORE_IPF
        WHERE SCORE_IPF.QVALUE < $IPF_FDR
        ) AS SCORE_IPF ON SCORE_IPF.FEATURE_ID = FEATURE.ID
    WHERE PRECURSOR.DECOY = 0 AND RUN.ID = $runID
    AND SCORE_IPF.PEPTIDE_ID+1 = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID
    ORDER BY transition_group_id, peak_group_rank;"
  } else if (runType == "DIA_Metabolomics") {
      query <- "SELECT PRECURSOR.ID AS transition_group_id,
         FEATURE.ID AS feature_id,
         FEATURE.EXP_RT AS RT,
         FEATURE_MS2.AREA_INTENSITY AS intensity,
         FEATURE.LEFT_WIDTH AS leftWidth,
         FEATURE.RIGHT_WIDTH AS rightWidth,
         SCORE_MS2.RANK AS peak_group_rank,
         SCORE_MS2.QVALUE AS m_score
         FROM PRECURSOR
         INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
         INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
         LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
         LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
         WHERE RUN.ID = $runID AND SCORE_MS2.QVALUE < $FDR AND PRECURSOR.DECOY = 0
         ORDER BY transition_group_id, peak_group_rank;"
  }
  query
}

#' Get precursor Info
#'
#' For each precursor in the table respective transition ids are fetched.
#' Order of transition is kept same as the order of their intensities in \code{\link{getTransitionsQuery}}.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-04-04
#' @param analytes (integer) A vector of integer that is searched in PRECURSOR.ID.
#' @param runType (char) This must be one of the strings "DIA_Proteomics", "DIA_IPF", "DIA_Metabolomics".
#' @return SQL query to be searched.
#' @seealso \code{\link{fetchPrecursorsInfo}}
#' @keywords internal
getPrecursorsQueryID <- function(analytes, runType = "DIA_Proteomics"){
  selectAnalytes <- paste0(" transition_group_id IN ('", paste(analytes, collapse="','"),"')")
  if (runType == "DIA_Proteomics" || runType == "DIA_IPF"){
      query <- paste0("SELECT PRECURSOR.ID AS transition_group_id,
                  TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
                  PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID AS peptide_id,
                  PEPTIDE.MODIFIED_SEQUENCE AS sequence,
                  PRECURSOR.CHARGE AS charge,
                  PRECURSOR.GROUP_LABEL AS group_label
                  FROM PRECURSOR
                  INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
                  INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
                  INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
                  WHERE ", selectAnalytes, "
                  ORDER BY peptide_id, transition_group_id, transition_id;")
  }
  else if (runType == "DIA_Metabolomics"){
      query <- paste0("SELECT PRECURSOR.ID AS transition_group_id,
                  TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id,
                  PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID AS peptide_id,
                  COMPOUND.SUM_FORMULA AS sum_formula,
                  COMPOUND.COMPOUND_NAME AS compound_name,
                  COMPOUND.ADDUCTS AS adducts,
                  PRECURSOR.CHARGE AS charge,
                  PRECURSOR.GROUP_LABEL AS group_label
                  FROM PRECURSOR
                  INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
                  INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID = PRECURSOR.ID
                  INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID
                  WHERE ", selectAnalytes, "
                  ORDER BY peptide_id, transition_group_id, transition_id;")
  }
  query
}


#' Get peptide scores
#'
#' For each peptide, its score, pvalue and qvalues are fetched across all runs.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-07-01
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_IPF", "DIA_Metabolomics".
#' @return SQL query to be searched.
#' @seealso \code{\link{getPeptideScores}}
#' @keywords internal
getPeptideQuery <- function(runType = "DIA_Proteomics"){
  if(runType == "DIA_Proteomics" || runType == "DIA_IPF"){
    query <- "SELECT DISTINCT SCORE_PEPTIDE.PEPTIDE_ID AS peptide_id,
  SCORE_PEPTIDE.RUN_ID AS run,
  SCORE_PEPTIDE.SCORE AS score,
  SCORE_PEPTIDE.PVALUE AS pvalue,
  SCORE_PEPTIDE.QVALUE AS qvalue
  FROM SCORE_PEPTIDE
  WHERE SCORE_PEPTIDE.CONTEXT = $CONTEXT;"
  } else if(runType == "DIA_Metabolomics"){
    query <- "SELECT DISTINCT COMPOUND.ID AS peptide_id,
  MIN(SCORE_MS2.QVALUE) AS minV,
  FEATURE.RUN_ID AS run,
  SCORE_MS2.SCORE AS score,
  SCORE_MS2.PVALUE AS pvalue,
  SCORE_MS2.QVALUE AS qvalue
  FROM SCORE_MS2
  INNER JOIN FEATURE ON FEATURE.ID = SCORE_MS2.FEATURE_ID
  INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID = FEATURE.PRECURSOR_ID
  INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID
  GROUP BY peptide_id, run;"
  }
  query
}

#' Get peptide scores
#'
#' For each peptide, its score, pvalue and qvalues are fetched across all runs.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-11-18
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_IPF", "DIA_Metabolomics".
#' @return SQL query to be searched.
#' @seealso \code{\link{getPeptideScores}}
#' @keywords internal
getPeptideQuery2 <- function(runType = "DIA_Proteomics"){
  query <- "SELECT DISTINCT SCORE_PEPTIDE.PEPTIDE_ID AS peptide_id,
  SCORE_PEPTIDE.RUN_ID AS run,
  SCORE_PEPTIDE.SCORE AS score,
  SCORE_PEPTIDE.PVALUE AS pvalue,
  SCORE_PEPTIDE.QVALUE AS qvalue
  FROM SCORE_PEPTIDE
  WHERE SCORE_PEPTIDE.CONTEXT = $CONTEXT AND SCORE_PEPTIDE.RUN_ID = $runID;"
  query
}

#' Get transitions from a SQLite file
#'
#' Query is generated to identify features and their transitions below a FDR cut-off from a run.
#' Order of transition intensity is kept same as the order of their Ids in \code{\link{getPrecursorsQuery}}.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-11-15
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_IPF", "DIA_Metabolomics".
#' @return SQL query to be searched.
#' @seealso \code{\link{fetchTransitionsFromRun}, \link{getPrecursorsQueryID}, \link{getPrecursorsQuery}}
#' @keywords internal
getTransitionsQuery <- function(runType = "DIA_Proteomics"){
  query <- "SELECT PRECURSOR.ID AS transition_group_id,
  FEATURE.ID AS feature_id,
  FEATURE.EXP_RT AS RT,
  FEATURE_TRANSITION.AREA_INTENSITY AS intensity,
  FEATURE.LEFT_WIDTH AS leftWidth,
  FEATURE.RIGHT_WIDTH AS rightWidth,
  SCORE_MS2.RANK AS peak_group_rank,
  SCORE_MS2.QVALUE AS m_score
  FROM PRECURSOR
  INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
  INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
  LEFT JOIN (
    SELECT FEATURE_ID, TRANSITION_ID, AREA_INTENSITY
    FROM FEATURE_TRANSITION
    ) AS FEATURE_TRANSITION ON FEATURE.ID = FEATURE_TRANSITION.FEATURE_ID
  INNER JOIN (
      SELECT FEATURE_ID, RANK, QVALUE
      FROM SCORE_MS2
      WHERE SCORE_MS2.QVALUE < $FDR
      ) AS SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
  WHERE PRECURSOR.DECOY = 0 AND RUN.ID = $runID
  ORDER BY transition_group_id, peak_group_rank, FEATURE_TRANSITION.TRANSITION_ID;"
  query
}


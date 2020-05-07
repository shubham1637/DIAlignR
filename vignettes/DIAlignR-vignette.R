## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installDIAlignR, eval=FALSE----------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("DIAlignR")

## ----loadDIAlignR-------------------------------------------------------------
library(DIAlignR)

## ----getDataPath--------------------------------------------------------------
dataPath <- system.file("extdata", package = "DIAlignR")

## ---- results=FALSE, message=FALSE, warning=FALSE-----------------------------
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
          "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
# For specific runs provide their names.
alignTargetedRuns(dataPath = dataPath, outFile = "test.csv", runs = runs, oswMerged = TRUE)
# For all the analytes in all runs, keep them as NULL.
alignTargetedRuns(dataPath = dataPath, outFile = "test.csv", runs = NULL, oswMerged = TRUE)

## ---- message=FALSE-----------------------------------------------------------
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
          "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
AlignObjLight <- getAlignObjs(analytes = 4618L, runs = runs, dataPath = dataPath, objType	= "light")
# First element contains names of runs, spectra files, chromatogram files and feature files.
AlignObjLight[[1]][, c("runName", "spectraFile")]
obj <- AlignObjLight[[2]][["4618"]][[1]][["AlignObj"]]
slotNames(obj)
names(as.list(obj))
AlignObjMedium <- getAlignObjs(analytes = 4618L, runs = runs, dataPath = dataPath, objType	= "medium")
obj <- AlignObjMedium[[2]][["4618"]][[1]][["AlignObj"]]
slotNames(obj)

## ---- fig.width=6, fig.align='center', fig.height=6, message=FALSE------------
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
 "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
AlignObj <- getAlignObjs(analytes = 4618L, runs = runs, dataPath = dataPath)
plotAlignedAnalytes(AlignObj, annotatePeak = TRUE)

## ---- fig.width=5, fig.align='center', fig.height=5, message=FALSE------------
library(lattice)
runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
 "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
AlignObjOutput <- getAlignObjs(analytes = 4618L, runs = runs,
                               dataPath = dataPath, objType = "medium")
plotAlignmentPath(AlignObjOutput)

## ---- fig.width=7, fig.align='center', fig.height=3.5, message=FALSE----------
data("XIC_QFNNTDIVLLEDFQK_3_DIAlignR")
XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]]
XICs.sm <- smoothXICs(XICs, type = "sgolay", samplingTime = 3.42, kernelLen = 9, polyOrd = 3)
plotXICgroup(XICs, Title = "Raw chromatograms")
plotXICgroup(XICs.sm, Title = "Smoothed chromatograms")

## -----------------------------------------------------------------------------
time <- lapply(XICs, `[[`, 1)
intensity <- lapply(XICs, `[[`, 2)
areaIntegrator(time, intensity, left = 5203.7, right = 5268.5, integrationType = "intensity_sum", baselineType = "base_to_base", fitEMG = FALSE)

## -----------------------------------------------------------------------------
sessionInfo()


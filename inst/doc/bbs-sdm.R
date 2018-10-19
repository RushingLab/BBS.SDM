## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  eval = FALSE,
  comment = "#>"
)

## ------------------------------------------------------------------------
#  library(BBS.tenstop)
#  library(BBS.SDM)

## ------------------------------------------------------------------------
#  bbs <- get_BBS10()

## ---- eval = FALSE-------------------------------------------------------
#  target_dir <- here::here("output") # change 'output' to the desired directory (relative to the current WD)

## ---- eval = FALSE-------------------------------------------------------
#  GetCorrData(bbs_raw = bbs, alpha = "FICR", path = target_dir,
#              start.year = 1972, end.year = 2014)

## ---- eval = FALSE-------------------------------------------------------
#  GetBioVars("FICR", path = target_dir)

## ---- eval = FALSE-------------------------------------------------------
#  GetInits(alpha = "FICR", path = target_dir)

## ---- eval = FALSE-------------------------------------------------------
#  RunMod(alpha = "FICR", path = target_dir, Parallel = TRUE, nI = 100, nB = 10, nC = 1)

## ---- eval = FALSE-------------------------------------------------------
#  ppc(alpha = "FICR", path = target_dir)

## ---- eval = FALSE-------------------------------------------------------
#  GetOccProb(alpha = "FICR", path = target_dir)
#  OccSummary(alpha = "FICR", path = target_dir)
#  GetIndices(alpha = "FICR", path = target_dir)


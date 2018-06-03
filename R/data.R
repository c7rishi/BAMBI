#' Backbone Dihedral Angles of Triose Phosphate Isomerase (8TIM)
#'
#' A dataset consisting of 490 pairs of backbone dihedral angles (in radian scale \eqn{[0, 2\pi)} )
#' \eqn{(\phi, \psi)} for the protein Triose Phosphate Isomerase (8TIM). The angles were obtained first by using
#' the DSSP software on the PDB file for 8TIM to get the backbone angles (in degrees),
#' and then by converting all angles into radians. Due to the presence of different secondary structures
#' (helices, sheets and loops) in the protein, the angular data show considerable variability, and is multimodal
#' with noticeably distinct clusters.
#'
#'
#' @format A data frame with 490 rows and 2 variables (backbone dihedral angles) phi and psi.
#' @source 8TIM PDB file: \url{http://www.rcsb.org/pdb/explore.do?structureId=8tim}.
#' @source DSSP software: \url{http://swift.cmbi.ru.nl/gv/dssp/}.
#'
#' @usage
#' data(tim8)

"tim8"



#' Saturna Island wind directions
#'
#' @description
#' A dataset consisting of 239 observations on wind direction in radians (original measurements were
#' in 10s of degrees), measured at Saturna Island, British Columbia,
#' Canada during October 1-10, 2016 (obtained from Environment Canada website). There was a severe storm
#' during October 4-7, which caused significant fluctuations among the wind directions. As a result the
#' angular data show a clear multimodality.
#'
#' @format A data frame with 239 rows and 2 columns; the column "angle" provides the angular direction (in radian)
#' and the column day provides the days on which the data points were collected (ranges between 1-10, corresponding to
#' October 1-10, 2016).
#' @source Environment Canada: \url{http://climate.weather.gc.ca/climate_data/data_quality_e.html}.
#' @source CBC news on the storm: \url{http://www.cbc.ca/news/canada/british-columbia/storm-bc-1.3795204}.
#'
#' @usage
#' data(wind)

"wind"

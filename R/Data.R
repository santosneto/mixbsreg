#' Antibody Level
#'
#' Data about a safety and immunogenecity study related to measles vaccines in Haiti.
#'
#' @docType data
#'
#' @usage data(haiti)
#'
#' @format A data frame with 330 observations on the following 4 variables.
#'
#'@details
#' \itemize{
#' \item IU neutralization antibody level.
#' \item EZ EZ is the type of vaccine used (0 if Schwartz and 1 if Edmonston-Zagreb).
#' \item HI HI is the level of the dosage (0 if medium and 1 if high).
#' \item FEM FEM is the gender (0 for male and 1 for female)
#'}
#'
#' @keywords datasets
#'
#' @references Moulton, L. H. and Halsey, N. A. (1995). A mixture model with detection limits for regression analyses of antibody response to vaccine. Biometrics, 51:1570–1578.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/8589241}{PubMed})
#'
#' @source \href{10.2307/2533289}{DOI}
#'
#' @examples
#' data(haiti)
#' hist(haiti$IU)
"haiti"

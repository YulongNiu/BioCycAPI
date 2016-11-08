##' BioCyc Database API - Get whole pathway list of a given species from BioCyc database.
##'
##' Get pathways from a given species. It may take more than 10 minutes to retrieve the xml file.
##' @title Get pathways
##' @inheritParams getCycTU
##' @return A vector of pathway ids.
##' @examples
##' ## Candida albicans SC5314 pathways
##' calTU <- getCycPathway('CALBI')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom xml2 read_xml xml_text xml_find_all
getCycPathway <- function(speID) {

  url <- paste0('http://biocyc.org/xmlquery?[x:x%3C-', speID, '^^pathways]')
  pathxml <- read_xml(url)

  pathVec <- xml_text(xml_find_all(pathxml, '//Pathway/@ID'))

  return(pathVec)
}

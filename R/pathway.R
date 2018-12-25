##' BioCyc Database API - Pathway
##'
##' getCycPathway(): Get whole pathways list,
##' getCycGenesfPathway(): Get genes from a given pathway ID.
##'
##' It may take more than 10 minutes to retrieve the xml file.
##' @title Pathway
##' @inheritParams getCycTU
##' @return
##'
##' getCycPathway(): A vector of pathway ids.
##' getCycGenesfPathway(): A vector of gene ids.
##'
##' @examples
##' ## Candida albicans SC5314 pathways
##' calTU <- getCycPathway('CALBI')
##'
##' ## genes of pathway "CALBI:PWY3B3-8"
##' genes <- getCycGenesfPathway('CALBI:PWY3B3-8')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml2 read_xml xml_text xml_find_all
##' @rdname pathway
##' @export
##'
getCycPathway <- function(speID) {

  url <- paste0('http://biocyc.org/xmlquery?[x:x%3C-', speID, '^^pathways]')
  pathxml <- read_xml(url)

  pathVec <- xml_text(xml_find_all(pathxml, '//Pathway/@ID'))

  return(pathVec)
}



##' @param pathID Pathway ID.
##' @importFrom xml2 read_xml xml_text xml_find_all
##' @rdname pathway
##' @export
getCycGenesfPathway <- function(pathID) {
  url <- paste0('http://websvc.biocyc.org/apixml?fn=genes-of-pathway&id=', pathID)
  genexml <- read_xml(url)

  geneVec <- xml_text(xml_find_all(genexml, '//Gene/@ID'))

  return(geneVec)
}


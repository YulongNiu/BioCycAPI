##' BioCyc Database API - Pathway
##'
##' \itemize{
##'   \item \code{getCycPathway()}: Get whole pathways list.
##'   \item \code{getCycGenesfPathway()}: Get genes from a given pathway ID.
##' }
##'
##' It may take more than 10 minutes to retrieve the xml file.
##' @title Pathway
##' @inheritParams getCycGenes
##' @return
##' \itemize{
##'   \item \code{getCycPathway()}:  A \code{tbl_df} contains 1st column is pathway ID and 2nd column is pathway annotation.
##'   \item \code{getCycGenesfPathway()}: A \code{list} indicates genes, proteins, and common names.
##' }
##'
##' @examples
##' ## Candida albicans SC5314 pathways
##' calTU <- getCycPathway('CALBI')
##'
##' ## genes of pathway "CALBI:PWY3B3-8"
##' genes <- getCycGenesfPathway('CALBI:PWY3B3-8')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom urltools url_encode
##' @importFrom xml2 read_xml xml_text xml_find_all
##' @importFrom magrittr %>%
##' @importFrom tibble tibble
##' @rdname pathway
##' @export
##'
getCycPathway <- function(speID) {

  pathxml <- paste0('http://biocyc.org/xmlquery?[x:x<-', speID, '^^pathways]') %>%
    url_encode %>%
    read_xml

  pathnodes <- pathxml %>%
    xml_find_all('/ptools-xml/Pathway')

  pathanno <- sapply(pathnodes, function(x){
    eachanno <- xml_find_all(x, 'common-name') %>%
      xml_text %>%
      ifelse(length(.) == 0, '', .)
    return(eachanno)
  })

  res <- tibble(pathID = xml_attr(pathnodes, 'ID'),
                pathAnno = pathanno)

  return(res)
}



##' @param pathID Pathway ID.
##' @importFrom xml2 read_xml xml_text xml_find_all
##' @rdname pathway
##' @export
getCycGenesfPathway <- function(pathID) {

  speID <- pathID %>%
    strsplit(split = ':', fixed = TRUE) %>%
    sapply('[[', 1)

  genexml <- paste0('http://websvc.biocyc.org/apixml?fn=genes-of-pathway&id=', pathID) %>%
    url_encode %>%
    read_xml

  res <- list(genes = genexml %>%
                xml_find_all('//Gene/@ID') %>%
                xml_text,
              proteins = genexml %>%
                xml_find_all('//Gene/product/Protein/@frameid') %>%
                xml_text %>%
                paste(speID, ., sep = ':'),
              names = genexml %>%
                xml_find_all('//Gene/common-name') %>%
                xml_text)

  return(res)
}


##' BioCyc Database Additional API - Get the NCBI taxonomy ID from a given BioCyc ID
##'
##' NCBI taxonomy ID is used as unique ID accoss cyc and BioCyc databases. This functions is used to get the corresponding NCBI Taxonomy ID from BioCyc. If BioCyc has no official NCBI taxonomy ID, it will return a character with length of 0 (it looks like "").
##' @title Get NCBI Taxonomy ID From BioCyc species ID
##' @inheritParams getCycGenes
##' @return The corresponding NCBI Taxonomy ID.
##' @examples
##' ## get human, mouse, and Ecoli NCBI taxonomy ID.
##' PhyloCyc2NCBI('HUMAN')
##' PhyloCyc2NCBI('MOUSE')
##' PhyloCyc2NCBI('ECOLI')
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml2 read_html xml_find_all xml_text
##' @importFrom magrittr %>%
##' @export
##'
PhyloCyc2NCBI <- function(speID) {
  ## USE: get the NCBI taxonomy ID from one BioCyc species ID
  ## INPUT: 'cycSpeID' is one BioCyc species ID
  ## OUTPUT: the NCBI taxnomy ID

  ## get cycID webpage
  web <- paste0('https://biocyc.org/', speID, '/organism-summary?object=', speID) %>%
    read_html

  ## get Taxonomy ID. The taxonomy ID is in the web-link like:
  ## <a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=237561" class="LINK">NCBI-Taxonomy:237561</a> </p> <P><P><CENTER><TABLE BORDER=1>
  taxID <- web %>%
    xml_find_all('//a[contains(@href, "www.ncbi.nlm.nih.gov/Taxonomy")]') %>%
    xml_text %>%
    strsplit(split = ':', fixed = TRUE) %>%
    sapply('[[', 2)

  return(taxID)
}

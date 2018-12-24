##' BioCyc Database Additional API - Get the NCBI taxonomy ID from a given BioCyc ID
##'
##' NCBI taxonomy ID is used as unique ID accoss cyc and BioCyc databases. This functions is used to get the corresponding NCBI Taxonomy ID from BioCyc. If BioCyc has no official NCBI taxonomy ID, it will return a character with length of 0 (it looks like "").
##' @title Get NCBI Taxonomy ID From BioCyc species ID
##' @param cycID The cyc species ID. The KEGG support multiple species ID, for example c('HUMAN', 'MOUSE', 'ECOLI').
##' @param n The number of CPUs or processors, and the default value is 4.
##' @return The corresponding NCBI Taxonomy ID in character vector.
##' @examples
##' ## get human, mouse, and Ecoli NCBI taxonomy ID with 2 threads
##' transPhyloCyc2NCBI(c('HUMAN', 'MOUSE', 'ECOLI'), n = 2)
##'
##' \dontrun{
##' ## transfer all BioCyc species ID to NCBI taxonomy ID
##' wBiocycSpe <- getCycPhylo(whole = TRUE)
##' wNCBISpe <- transPhyloCyc2NCBI(wBiocycSpe[, 1])
##' }
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom RCurl getURL
##' @importFrom doMC registerDoMC
##' @importFrom foreach foreach %dopar%
##' @importFrom KEGGAPI getcontent
##' @export
##'
transPhyloCyc2NCBI <- function(cycID, n = 4){

  registerDoMC(n)

  getSingleTax <- function(cycSpeID) {
    ## USE: get the NCBI taxonomy ID from one BioCyc species ID
    ## INPUT: 'cycSpeID' is one BioCyc species ID
    ## OUTPUT: the NCBI taxnomy ID

    ## get cycID webpage
    cycLink <- paste('http://biocyc.org/', cycSpeID, '/organism-summary?object=', cycSpeID, sep = '')
    cycWeb <- getURL(cycLink)

    ## get Taxonomy ID. The taxonomy ID is in the web-link like 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=593907'
    taxIDLink <- gregexpr('wwwtax\\.cgi\\?mode=Info\\&id=\\d+', cycWeb)
    taxIDLink <- getcontent(cycWeb, taxIDLink[[1]])
    taxID <- gregexpr('\\d+', taxIDLink)
    taxID <- getcontent(taxIDLink, taxID[[1]])

    return(taxID)
  }

  cycTax <- foreach(i = 1:length(cycID), .combine = c) %dopar% {
   print(paste0('It is running ', i, ' in a total number of ', length(cycID), '.'))
    singleTaxID <- getSingleTax(cycID[i])
    names(singleTaxID) <- cycID[i]
    return(singleTaxID)
  }

  return(cycTax)

}

##' BioCyc Database API - Get transcription unit (TU) from gene
##'
##' Get TU from a given BioCyc gene ID.
##' If the given gene has no TU, "NULL" will be returned. If the "evidence" is set to TRUE, a list will return.
##' @title Get TU from ID.
##' @param geneID A BioCyc gene.
##' @param evidence  Logical value indicates whether to return the evidence value.
##' @return A \code{list} contains TUs or NULL
##' @examples getCycTUfGene('ECOLI:EG10102')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml read_xml
##' @export
##'
getCycTUfGene <- function(geneID, evidence = FALSE) {

  ## read in TU information XML
  url <- paste0('http://biocyc.org/apixml?fn=transcription-units-of-gene&id=', geneID, '&detail=low') %>%
    url_encode

  tuxml <- read_xml(url)

  res <- list(TU = tuxml %>%
                xml_find_all('/ptools-xml/Transcription-Unit/@ID') %>%
                xml_text,
              evidence = tuxml %>%
                xml_find_all('/ptools-xml/Transcription-Unit/evidence/Evidence-Code/@ID') %>%
                xml_text)

  return(res)
}


##' BioCyc Database API - Get TU information from BioCyc database.
##'
##' Get TU information including genes in TU and evidence.
##' There is another way to get the genes "http://biocyc.org/apixml?fn=transcription-unit-genes&id=ECOLI:TU0-42328&detail=full".
##' @title Get one TU information
##' @param TUID A BioCyc TU ID with the length of 1.
##' @param speID The BioCyc species ID, for example "ECOLI" is for "Escherichia coli K-12 substr. MG1655".
##' @return A list of TU information.
##' @examples
##' getCycTUInfo('ECOLI:TU0-6636')
##' getCycTUInfo('ECOLI:TU00260')
##' getCycTUInfo('SMUT210007:TU1FZX-806')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom urltools url_encode
##' @importFrom xml read_xml xml_text xml_find_all
##' @export
##'
getCycTUInfo <- function(TUID) {

  ## species identity
  speid <- TUID %>%
    strsplit(split = ':', fixed = TRUE) %>%
    sapply('[[', 1)

  ## read in gene information XML
  url <- paste0('http://biocyc.org/apixml?fn=transcription-units-of-gene&id=', TUID, '&detail=full') %>%
    url_encode

  tuxml <- read_xml(url)

  ## genes
  genes <- tuxml %>%
    xml_find_all('//component/Gene/@frameid') %>%
    xml_text %>%
    paste(speid, ., sep = ':')

  ## terminator
  terminator <- list(id = tuxml %>%
                       xml_find_all('//Terminator/@ID') %>%
                       xml_text,
                     right = tuxml %>%
                       xml_find_all('//Terminator/right-end-position') %>%
                       xml_text,
                     left =  tuxml %>%
                       xml_find_all('//Terminator/left-end-position') %>%
                       xml_text,
                     evidence = tuxml %>%
                       xml_find_all('//Terminator/evidence/Evidence-Code/@ID') %>%
                       xml_text)

  ## promoter
  promoter <- list(id = tuxml %>%
                     xml_find_all('//Promoter/@ID') %>%
                     xml_text,
                   right = tuxml %>%
                     xml_find_all('//Promoter/right-end-position') %>%
                     xml_text,
                   left =  tuxml %>%
                     xml_find_all('//Promoter/left-end-position') %>%
                     xml_text,
                   evidence = tuxml %>%
                     xml_find_all('//Promoter/evidence/Evidence-Code/@ID') %>%
                     xml_text)


  ## TU evidence
  TUEv <- tuxml %>%
    xml_find_all('//Transcription-Unit/evidence/Evidence-Code/@ID') %>%
    xml_text

  ## merge
  res <- list(genes = genes,
              promoter = promoter,
              terminator = terminator,
              evidence = TUEv,
              url = url)

  return(res)
}


##' BioCyc Database API - Get whole transcription unit list of a given species from BioCyc database.
##'
##' Get transcription units from a given species. It may take more than 10 minutes to retrieve the xml file.
##'
##' @title TU
##' @param speID The BioCyc species ID, for example "ECOLI" is for "Escherichia coli K-12 substr. MG1655".
##' @return A vector of TU ids.
##' @examples
##' ## get Streptococcus mutans UA159 TU
##' smTU <- getCycTU('SMUT210007')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml2 read_xml xml_text xml_find_all
##' @importFrom urltools url_encode
##' @importFrom magrittr %>%
##' @export
##'
getCycTU <- function(speID){

  tuxml <- paste0('http://websvc.biocyc.org/xmlquery?[x:x<-', speID, '^^Transcription-Units]') %>%
    url_encode %>%
    read_xml

  res <- tuxml %>%
    xml_find_all('//Transcription-Unit/@ID') %>%
    xml_text

  return(res)
}

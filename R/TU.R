##' BioCyc Database API - Get transcription unit (TU) from gene
##'
##' Get TU from a given BioCyc gene ID.
##' If the given gene has no TU, "NULL" will be returned. If the "evidence" is set to TRUE, a list will return.
##' @title Get TU from ID.
##' @param geneID A BioCyc gene.
##' @param speID The BioCyc species ID, for example "ECOLI" is for "Escherichia coli K-12 substr. MG1655".
##' @param evidence  Logical value indicates whether to return the evidence value.
##' @return A vector contains TUs or NULL
##' @examples getCycTUfGene('EG10102', 'ECOLI')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse
##' @export
##'
##' 
getCycTUfGene <- function(geneID, speID, evidence = FALSE) {

  ## read in TU information XML
  url <- paste0('http://biocyc.org/apixml?fn=transcription-units-of-gene&id=', speID, ':', geneID, '&detail=low')
  TUInfoXML <- xmlRoot(xmlTreeParse(url))

  ## get TU
  TUID <- xmlNodeAttr(TUInfoXML, '/ptools-xml/Transcription-Unit', 'frameid')
  TUID <- testLen(TUID, NULL, TUID)

  return(TUID)
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
##' getCycTUInfo('TU0-6636', 'ECOLI')
##' getCycTUInfo('TU00260', 'ECOLI')
##' getCycTUInfo('TUC7Z-43', 'SMUT210007')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom XML xmlRoot xmlTreeParse
##' @export
##'
getCycTUInfo <- function(TUID, speID) {

  ## read in gene information XML
  url <- paste('http://biocyc.org/apixml?fn=transcription-units-of-gene&id=', speID, ':', TUID, '&detail=full', sep = '')
  TUInfoXML <- xmlRoot(xmlTreeParse(url))

  ## genes in the TU
  geneIDs <- xmlNodeAttr(TUInfoXML, '//component/Gene', 'frameid')

  ## terminator
  terminatorName <- xmlNodeAttr(TUInfoXML, '//component/Terminator', 'frameid')
  terminatorName <- testLen(terminatorName, NULL, terminatorName)
  ## right end
  rightPos <- xmlNodeVal(TUInfoXML, '//component/Terminator/right-end-position')
  rightPos <- testLen(rightPos, NULL, rightPos)
  ## left end
  leftPos <- xmlNodeVal(TUInfoXML, '//component/Terminator/left-end-position')
  leftPos <- testLen(leftPos, NULL, leftPos)
  ## teminator evidence
  TMEv <- xmlNodeAttr(TUInfoXML, '//component/Terminator/evidence/Evidence-Code', 'frameid')
  TMEv <- testLen(TMEv, NULL, TMEv)
  terminator <- list(name = terminatorName,
                     rightPos = rightPos,
                     leftPos = leftPos,
                     Ev = TMEv)

  ## promoter
  promoterName <- xmlNodeAttr(TUInfoXML, '//component/Promoter', 'frameid')
  promoterName <- testLen(promoterName, NULL, promoterName)
  ## promoter common name
  promoterComName <- xmlNodeVal(TUInfoXML, '//component/Promoter/common-name')
  promoterComName <- testLen(promoterComName, NULL, promoterComName)
  ## right end
  rightPos <- xmlNodeVal(TUInfoXML, '//component/Promoter/right-end-position')
  rightPos <- testLen(rightPos, NULL, rightPos)
  ## left end
  leftPos <- xmlNodeVal(TUInfoXML, '//component/Promoter/left-end-position')
  leftPos <- testLen(leftPos, NULL, leftPos)
  ## promoter evidence
  PMEv <- xmlNodeAttr(TUInfoXML, '//component/Promoter/evidence/Evidence-Code', 'frameid')
  PMEv <- testLen(PMEv, NULL, PMEv)
  promoter <- list(name = promoterName,
                   comName = promoterComName,
                   rightPos = rightPos,
                   leftPos = leftPos,
                   Ev = PMEv)
  
  ## TU evidence
  TUEv <- xmlNodeAttr(TUInfoXML, '//Transcription-Unit/evidence/Evidence-Code', 'frameid')
  TUEv <- testLen(TUEv, NULL, TUEv)

  ## merge
  cycTU <- list(geneIDs = geneIDs,
                terminator = terminator,
                promoter = promoter,
                TUEv = TUEv,
                url = url)

  return(cycTU)

}


##' BioCyc Database API - Get whole transcription unit list of a given species from BioCyc database.
##' Get transcription units from a given species. It may take more than 10 minutes to retrieve the xml file.
##'
##' @title TU
##' @param speID The BioCyc species ID, for example "ECOLI" is for "Escherichia coli K-12 substr. MG1655".
##' @return A vector of TU ids.
##' @examples
##' ## get Streptococcus mutans UA159 TU
##' smTU <- getCycTU('SMUT210007')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom xml2 read_xml xml_text xml_find_all
##' @export
getCycTU <- function(speID){

  url <- paste0('http://websvc.biocyc.org/xmlquery?[x:x%3C-', speID, '^^Transcription-Units]')
  TUxml <- read_xml(url)

  TUVec <- xml_text(xml_find_all(TUxml, '//Transcription-Unit/@ID'))

  return(TUVec)
}

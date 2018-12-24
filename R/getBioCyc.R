##' BioCyc Database API - Get species information from BioCyc.
##'
##' Get the BioCyc species information including the BioCyc species ID and the Latin name.
##' @title Get species from BioCyc
##' @param speList The species list that is a vector like 'c("HUMAN", "ECOLI", "ZMOB579138")'. The input speList should be consistent with the parameter 'speType'.
##' @param speType It supports two types: "BioCyc" and "regexpr".
##' BioCyc type is the BioCyc species ID, for exmaple "HUMAN" is the BioCyc ID for the Homo sapiens.The "regexpr" is used for regulare expression search with the Latin name for example "Escherichia coli".
##' @param whole Whether or not get the whole BioCyc species list,
##' and the default value is FALSE.
##' @return Matrix of species information.
##' @examples
##' ## search species list from BioCyc ID
##' getCycPhylo(c('HUMAN', 'ECOLI', 'ZMOB579138'), speType = 'BioCyc')
##'
##' ## search species whose names include 'Escherichia coli'
##' getCycPhylo('Escherichia coli', speType = 'regexpr')
##'
##' \dontrun{
##' ## get whole BioCyc species information table
##' getCycPhylo(whole = TRUE)}
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml2 read_xml xml_find_all xml_attrs xml_text xml_children
##' @export
##'
getCycPhylo <- function(speList, speType = 'BioCyc', whole = FALSE) {

  if (!(speType %in% c('BioCyc', 'regexpr'))) {
    stop('"speType" now only supports "BioCyc" and "regexpr".')
  } else {}

  ## read in the whole biocyc XML file
  cycSpeXML <- read_xml('http://biocyc.org/xmlquery?dbs')

  ## get each species
  cycSpe <- xml_find_all(cycSpeXML, '//PGDB')

  ## get biocyc ID for each species
  cycSpeID <- t(sapply(cycSpe, xml_attrs))

  ## get latin name of each species
  cycLatin <- lapply(cycSpe, function(x) {
    eachLatin <- xml_text(xml_children(x))
    return(eachLatin)
  })

  cycLatin <- sapply(cycLatin, paste, collapse = ' ')
  cycSpeMat <- cbind(cycSpeID[, 1], cycLatin, cycSpeID[, 2])
  colnames(cycSpeMat) <- c('BioCycID', 'LatinName', 'Version')

  if(!whole) {
    if (speType == 'BioCyc') {
      cycSpeMat <- cycSpeMat[cycSpeMat[, 1] %in% speList, , drop = FALSE]
    }
    else if (speType == 'regexpr') {
      cycSpeMat <- cycSpeMat[grep(speList, cycSpeMat[, 2]), , drop = FALSE]
    }
  } else {}

  return(cycSpeMat)
}


##' BioCyc Database API - Get the whole genes/proteins list from a given species
##'
##' Get the BioCyc gene ID or protein ID list of a given species.
##' @title Get whole genes/proteins list
##' @param speID The BioCyc species ID, for example "ECOLI" is for "Escherichia coli K-12 substr. MG1655".
##' @param type Get the "genes" or "proteins", and the default value is "genes".
##' @return A vector of genes of proteins with BioCyc ID.
##' @examples
##' getCycGenes('ECOLI')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml2 read_xml
##' @importFrom magrittr %>%
##' @export
##'
##'
getCycGenes <- function(speID, type = 'genes'){

  ## read in the whole biocyc XML file
  cycxml <- paste0('http://biocyc.org/xmlquery?query=[x:x<-', speID, '^^', type, ']&detail=none') %>%
    URLencode %>%
    read_xml

  ## select proteins or genes
  ## <Gene resource="getxml?ECOLI:G7553" orgid="ECOLI" frameid="G7553"/>
  ## remove <Gene resource="getxml?ECOLI:BC-2.2.2" orgid="ECOLI" frameid="BC-2.2.2" class="true"/>
  nd <- ifelse(type == 'genes', 'Gene', 'Protein')

  res <- nd %>%
    paste0(., '[not(@class)]') %>%
    xml_find_all(cycxml, .) %>%
    xml_attr(., 'frameid') %>%
    paste(speID, ., sep = ':')

  return(res)
}

##' BioCyc Database API - Get gene information from BioCyc database.
##'
##' The gene information from BioCyc including genome location and gene name information. Some genes in BioCyc may do not have common names or accession names In this circumstance, "NULL" will be return
##' @title Get one gene information from BioCyc.
##' @param geneID A BioCyc gene ID with the length of 1.
##' @param speID The BioCyc species ID, for example "ECOLI" is for "Escherichia coli K-12 substr. MG1655".
##' @return A list.
##' @examples
##' ## get "atpE" gene information from Ecoli K-12 MG1655 strain.
##' getCycGeneInfo('EG10102', 'ECOLI')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom XML xmlRoot xmlTreeParse getNodeSet xmlName
##' @export
##'
getCycGeneInfo <- function(geneID, speID){

  ## read in gene information XML
  url <- paste('http://biocyc.org/getxml?', speID, ':',geneID, '&detail=full', sep = '')
  geneInfoXML <- xmlRoot(xmlTreeParse(url))

  ## location in genome
  ## direction
  dirTrans <- xmlNodeVal(geneInfoXML, '//transcription-direction')
  ## right end
  rightPos <- xmlNodeVal(geneInfoXML, '//right-end-position')
  ## left end
  leftPos <- xmlNodeVal(geneInfoXML, '//left-end-position')
  locTrans <- list(dirTrans = dirTrans,
                   rightPos = rightPos,
                   leftPos = leftPos)

  ## gene names, some genes may not have names
  ## common name
  comName <- xmlNodeVal(geneInfoXML, '//common-name')
  comName <- testLen(comName, NULL, comName)

  ## accession, the accession name of node is like 'accession-1', 'accession-2'
  allNodeName <- sapply(getNodeSet(geneInfoXML, '//*'), xmlName)
  accNodeName <- allNodeName[grepl('^accession-\\d+', allNodeName)]
  ## some genes may have no accession name
  if (length(accNodeName) == 0) {
    accName = NULL
  } else {
    accNodePath <- paste('//', accNodeName, sep = '')
    accName <- unname(sapply(accNodePath, xmlNodeVal, xmlFile = geneInfoXML))
  }
  geneName <- list(comName = comName,
                   accName = accName)

  ## merge
  cycGene <- list(locTrans = locTrans,
                  geneName = geneName,
                  url = url)

  return(cycGene)

}


##' Translate KEGG ID to BioCyc ID.
##'
##' Translate the KEGG gene ID to BioCyc gene ID is tricky. In BioCyc, the gene names is not in a uniform; some of them use symbol like "dnaK", but some use the KEGGID. For genes, if symbol is given, we use 'http://biocyc.org/xmlquery?query=[x:x<-ECOLI^^genes,x^name="atpA"]&detail=full'. There are two circumstances that will return "0": one is that  BioCyc database may marker some genes as "Pseudo-Genes", and the other is different gene symbols in KEGG and BioCyc. For proteins, we at first transfer KEGG gene IDs to UniProt IDs, and then to BioCyc gene IDs.
##' 
##' @title Transfer KEGG ID to BioCyc ID.
##' @param KEGGID Only one KEGG ID
##' @param speKEGGID Species BioCyc ID.
##' @param speCycID Species KEGG ID
##' @param type 'gene' or 'protein'
##' @return The BioCyc gene ID or "0", if gene is not found.
##' @examples
##' ## symbol is "atpD"
##' transGeneIDKEGG2Cyc('b3732', 'eco', 'ECOLI')
##' 
##' ## symbol is "SMU_408" but the first annotation word is "permease"
##' transGeneIDKEGG2Cyc('SMU_408', 'smu', 'SMUT210007')
##' 
##' ## It will return "0" because of the symbol 'atpE_H' from KEGG.
##' ## The symbol in BioCyc is 'atpE/H'.
##' transGeneIDKEGG2Cyc('Bd0010', 'bba', 'BBAC264462')
##' 
##' ## retrieve protein
##' transGeneIDKEGG2Cyc('b0001', 'eco', 'ECOLI', type = 'protein')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml2 read_xml xml_children xml_attr
##' @importFrom KEGGAPI webTable convKEGG
##' @export
##'
##' 
transGeneIDKEGG2Cyc <- function(KEGGID, speKEGGID, speCycID, type = 'gene') {


  KEGGID2symbol <- function(KEGGIDtry, speKEGGIDtry) {
    ## transfer KEGG ID to symbol, if it has one; otherwise, we just use the KEGGID
    ## To save the gene symbol at most, the return gene information is split by ";" and ",". All elements, expect one contains space, are returned.
    ## USE: try to convert KEGG gene ID to symbole from KEGG gene annoation. If the first 'proper word' is all in lower case, return it plus KEGGID. If no 'proper word', KEGGID returns.
    ## INPUT: 'KEGGIDtry' is the KEGG gene ID. 'speKEGGIDtry' is the KEGG species ID.
    ## OUTPUT: A vector (length may be bigger than 1)
    ## EXAMPLE: KEGGID2symbol('SMU_23', 'smu')
    ## EXAMPLE: KEGGID2symbol('SMU_t02', 'smu')
    ## EXAMPLE: KEGGID2symbol('SMU_24', 'smu')
    ## EXAMPLE: KEGGID2symbol('SMU_20', 'smu')
    ## EXAMPLE: KEGGID2symbol('SMU_408', 'smu')
    KEGGsymTable <- webTable(paste('http://rest.kegg.jp/list/', speKEGGIDtry, ':', KEGGIDtry, sep = ''), ncol = 2)
    KEGGsym <- KEGGsymTable[1, 2]

    ## split by '; ' and ', '
    ## if KEGGsym is '', then NULL will return. c('test', NULL) is equal to 'test'
    KEGGEle <- unlist(strsplit(KEGGsym, split = '; ', fixed = TRUE))
    KEGGEle <- unlist(strsplit(KEGGEle, split = ', ', fixed = TRUE))

    ## remove elements with space
    KEGGEle <- KEGGEle[!grepl(' ', KEGGEle)]

    ## combine KEGGIDtry
    KEGGEle <- c(KEGGIDtry, KEGGEle)

    return(KEGGEle)
  }

  TrySymConv <- function(symboltry, speCycIDtry) {
    ## USE: try to convert KEGG symbol to Biocyc gene ID
    ## INPUT: 'symboltry' is the gene symbol extract from KEGG. 'speCycIDtry' is Biocyc species ID.
    ## OUTPU: converted BioCyc gene ID. '"0"' will return, if not found.
    ## EXAMPLE: TrySymConv('dnaA', 'SMUT210007')
    ## EXAMPLE: TrySymConv('SMU_06', 'SMUT210007')

    ## ## old code for xml, but will not work because "%" --> "%25"
    ## url <- paste0('http://biocyc.org/xmlquery?query=[x:x<-', speCycIDtry, '^^genes,x^name','="', symboltry, '"]&detail=full')
    ## geneXML <- xmlRoot(xmlTreeParse(url))
    ## cycID <- xmlNodeAttr(geneXML, '/ptools-xml/Gene', 'frameid')
    ## cycID <- testLen(cycID, '0', cycID)

    ## try symbol
    url <- paste0('http://websvc.biocyc.org/xmlquery?query=[x:x%3C-',speCycIDtry, '^^genes,x^name%3D%22', symboltry, '%22]&detail=full')
    geneXML <- read_xml(url)
    geneXMLChild <- xml_children(geneXML)

    ## check if xml returns geneID
    if (length(geneXMLChild) < 2) {
      cycID = '0'
    } else {
      geneXMLChildCont <- geneXMLChild[[2]]
      cycID <- xml_attr(geneXMLChildCont, 'frameid', default = '0')
    }

    cycIDList <- list(cycID = cycID, url = url)
    ## # also get protein
    ## cycID <- xmlNodeAttr(geneXML, '//Protein', 'frameid')

    return(cycIDList)
  }

  if (type == 'gene') {
    ## convert to symbol
    symbol <- KEGGID2symbol(KEGGID, speKEGGID)
    for (i in 1:length(symbol)) {
      cycIDList <- TrySymConv(symbol[i], speCycID)
      if (cycIDList$cycID != '0') {
        break
      } else {}
    }
  }
  else if (type == 'protein') {
    # transfer KEGG ID to unipro ID
    standKEGGID <- paste(speKEGGID, KEGGID, sep = ':')
    uniproID <- convKEGG('uniprot', standKEGGID, convertType = 'identity', n = 1)
    uniproID <- uniproID[1, 2]

    # 'up:Q8DWN9' --> 'Q8DWN9'
    uniproID <- sapply(strsplit(uniproID, split = ':', fixed = TRUE), '[[', 2)
    url <- paste0('http://websvc.biocyc.org/', speCycID, '/foreignid?ids=Uniprot:', uniproID)
    cycIDMat <- webTable(url, ncol = 1)
    if (cycIDMat[2, 1] == 1) {
      cycID <- cycIDMat[3, 1]
    }
    else if (cycIDMat[2, 1] == 0) {
      cycID <- NULL
    }

    cycID <- testLen(cycID, '0', cycID)
    cycIDList <- list(cycID = cycID, url = url)
  }

  return(cycIDList)
}


## Transfer KEGG ID to BioCyc ID.
##
## Tranfser KEGG ID to BioCyc ID by the route, whole species gene list --> unit gene information --> match KEGG ID
## @title Transfer genes ID from KEGG to BioCyc
## @param KEGGID A vector of KEGG gene ID.
## @param speID Species BioCyc ID.
## @param n The number of CPUs or processors, and the default value is 4.
## @param type Only support genes
##  @return A vector of BioCycID
## # EG10098 EG10101
## # "b3734" "b3732"
## # @examples transGeneIDKEGG2Cyc(c('b3734', 'b3732'), 'ECOLI')
## # @author Yulong Niu \email{yulong.niu@@hotmail.com}
##  @importFrom doMC registerDoMC
##  @importFrom foreach foreach
##  @export
##
## transGeneIDKEGG2Cyc <- function(KEGGID, speID, n = 4, type = 'genes'){

##   require(foreach)
##   require(doMC)

##   registerDoMC(n)

##   # whole list of given species
##   geneList <- getCycGenesList(speID)

##   # whole accession number
##   accList <- foreach(i = 1:length(geneList)) %dopar% {
##     print(paste('It is running ', i, '.', sep = ''))
##     geneInfo <- getCycGeneInfo(geneList[i], speID)
##     accInfo <- geneInfo$geneName$accName
##     # some genes have no accession names
##     if (!is.null(accInfo)) {
##       names(accInfo) <- rep(geneList[i], length(accInfo))
##     } else {}
##     return(accInfo)
##   }

##   accList <- unlist(accList)

##   KEGGList <- accList[accList %in% KEGGID]

##   return(KEGGList)


## }


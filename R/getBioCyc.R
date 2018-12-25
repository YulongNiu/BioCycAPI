##' BioCyc Database API - Get species information from BioCyc.
##'
##' Get the BioCyc species information including the BioCyc species ID and the Latin name.
##' @title Get species from BioCyc
##' @param speList A \code{character vector} contains species list like 'c("HUMAN", "ECOLI", "ZMOB579138")'. It should be consistent with the parameter `speType`.
##' @param speType A \code{string} "BioCyc" or "regexpr". "BioCyc" means it is the BioCyc species ID, for example "HUMAN" is the BioCyc ID for the Homo sapiens.The "regexpr" is used for regulare expression search with the Latin name for example "Escherichia coli".
##' @param whole A \code{logic} whether or not get the whole BioCyc species list, and the default value is FALSE.
##' @param n A \code{integer} refers to the number of cores.
##' @return A \code{tbl_df} object contains BioCycID, Latin name, and version.
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
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom tibble tibble
##' @importFrom xml2 read_xml xml_find_all xml_attrs xml_text
##' @importFrom magrittr %>% %<>%
##' @importFrom dplyr filter
##' @export
##'
getCycPhylo <- function(speList, speType = 'BioCyc', whole = FALSE, n = 2) {

  if (!(speType %in% c('BioCyc', 'regexpr'))) {
    stop('"speType" now only supports "BioCyc" and "regexpr".')
  } else {}

  ## read in the whole biocyc XML file
  cycxml <- read_xml('http://biocyc.org/xmlquery?dbs')

  ## get each species
  cycspe <- cycxml %>%
    xml_find_all('//PGDB')

  ##~~~~~~~~~~~find species Latin name~~~~~~~~~~~~~~~~
  registerDoParallel(cores = n)
  latin <- foreach(i = seq_along(cycspe), .combine = c) %dopar% {
    eachname <- cycspe[i] %>%
      xml_find_all('.//text()') %>%
      xml_text %>%
      paste(collapse = ' ')
  }
  ## stop multiple cores
  stopImplicitCluster()
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  res <- tibble(BioCycID = xml_attr(cycspe, 'orgid'),
                LatinName = latin,
                Version = xml_attr(cycspe, 'version'))

  ## select by ID or regexpr
  if(!whole) {
    if (speType == 'BioCyc') {
      res %<>% filter(BioCycID %in% speList)
    }
    else if (speType == 'regexpr') {
      res %<>% filter(grepl(speList, LatinName))
    }
  } else {}

  return(res)
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
##' @importFrom xml2 read_xml xml_find_all xml_attr
##' @importFrom magrittr %>%
##' @importFrom urltools url_encode
##' @export
##'
##'
getCycGenes <- function(speID, type = 'genes'){

  ## read in the whole biocyc XML file
  cycxml <- paste0('http://biocyc.org/xmlquery?query=[x:x<-', speID, '^^', type, ']&detail=none') %>%
    url_encode %>%
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
##' @return A list.
##' @examples
##' ## get "atpE" gene information from Ecoli K-12 MG1655 strain.
##' getCycGeneInfo('ECOLI:EG10102')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml2 read_xml xml_find_all xml_text
##' @importFrom magrittr %>%
##' @importFrom urltools url_encode
##' @export
##'
getCycGeneInfo <- function(geneID){

  speID <- geneID %>%
    strsplit(split = ':', fixed = TRUE) %>%
    sapply('[[', 1)

  ## read in gene information XML
  url <- paste0('http://biocyc.org/getxml?', geneID, '&detail=full') %>%
    url_encode
  genexml <- read_xml(url)

  ##~~~~~~~~~~~~~~~~~~~~location in genome~~~~~~~~~~~~~~~~~~~~~~~
  loc <- list(direction = genexml %>%
                xml_find_all('//transcription-direction') %>%
                xml_text,
              right = genexml %>%
                xml_find_all('//right-end-position') %>%
                xml_text,
              left = genexml %>%
                xml_find_all('//left-end-position') %>%
                xml_text,
              centisome = genexml %>%
                xml_find_all('//centisome-position') %>%
                xml_text)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##~~~~~~~~~~~~~~~~gene names, some genes may not have names~~~~~
  ## select contains 'accession*', 'synonym', 'common-name'
  name <- genexml %>%
    xml_find_all('//*[contains(name(), "accession")] | //synonym | //common-name') %>%
    xml_text

  protein <- genexml %>%
    xml_find_all('//Protein/@frameid') %>%
    xml_text %>%
    paste(speID, ., sep = ':')
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##~~~~~~~~~~~~~~~~TU~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TU <- genexml %>%
    xml_find_all('//Transcription-Unit/@frameid') %>%
    xml_text %>%
    paste(speID, ., sep = ':')
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## merge
  res <- list(loc = loc,
              name = name,
              protein = protein,
              TU = TU,
              url = url)

  return(res)
}


##' Convert KEGG ID to BioCyc ID.
##'
##' Convert the KEGG gene ID to BioCyc gene ID is tricky. In BioCyc, the gene name is not in a uniform; some of them use symbol like "dnaK", but some use the KEGG ID. For genes, if symbol is given, we use 'http://biocyc.org/xmlquery?query=[x:x<-ECOLI^^genes,x^name="atpA"]&detail=full'. There are two circumstances that will return "0": one is that  BioCyc database may marker some genes as "Pseudo-Genes", and the other is different gene symbols in KEGG and BioCyc. For proteins, we at first transfer KEGG gene IDs to UniProt IDs, and then to BioCyc gene IDs.
##'
##' @title Transfer KEGG ID to BioCyc ID.
##' @param geneKEGG The gene symbol extract from KEGG.
##' @param speBioCyc Species KEGG ID
##' @param type 'gene' or 'protein'
##' @return The BioCyc gene ID or "0", if gene is not found.
##' @examples
##' ## symbol is "atpD"
##' ConvGeneKEGG2BioCyc('eco:b3732', 'ECOLI')
##'
##' ## symbol is "SMU_408" but the first annotation word is "permease"
##' ConvGeneKEGG2BioCyc('smu:SMU_408', 'SMUT210007')
##'
##' ## It will return "0" because of the symbol 'atpE_H' from KEGG.
##' ## The symbol in BioCyc is 'atpE/H'.
##' ConvGeneKEGG2BioCyc('bba:Bd0010', 'BBAC264462')
##'
##' ## retrieve protein
##' ConvGeneKEGG2BioCyc('eco:b0001', 'ECOLI', type = 'protein')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml2 read_xml xml_children xml_attr
##' @importFrom KEGGAPI webTable convKEGG
##' @importFrom magrittr %>%
##' @export
##'
ConvGeneKEGG2BioCyc <- function(geneKEGG, speBioCyc, type = 'gene') {

  if (type == 'gene') {
    ## convert to symbol
    symbol <- KEGG2Sym(geneKEGG)
    for (i in seq_along(symbol)) {
      res <- Sym2BioCyc(symbol[i], speBioCyc)
      if (length(res$gene) != 0) {
        break
      } else {}
    }
  }
  else if (type == 'protein') {
    ## strategy
    ## transfer KEGG ID to unipro ID
    uniproID <- geneKEGG %>%
      convKEGG('uniprot', ., convertType = 'identity', n = 1) %>%
      .[1, 2] %>%
      strsplit(split = ':', fixed = TRUE) %>%
      sapply('[[', 2) ## 'up:Q8DWN9' --> 'Q8DWN9'

    ## get BioCyc ID
    url <- paste0('https://websvc.biocyc.org/', speBioCyc, '/foreignid?ids=Uniprot:', uniproID)

    cycID <- url %>%
      webTable(ncol = 1) %>%
      as.character %>%
      .[length(.)]

    res <- list(protein = ifelse(cycID == '0', character(0), cycID),
                url = url)
  }

  return(res)
}


##' KEGG gene ID to symbol
##'
##' Try to transfer KEGG ID to symbol, if it has one; otherwise, we just use the `KEGGID`. The return gene information is split by ";" and ",". All elements, expect ones contain space, are returned.
##'
##' @title KEGG gene ID to symbol.
##' @inheritParams ConvGeneKEGG2BioCyc
##' @return A \code{character vector} including the raw KEGGID.
##' @examples
##' KEGG2Sym('smu:SMU_23')
##'
##' KEGG2Sym('smu:SMU_t02')
##'
##' KEGG2Sym('smu:SMU_24')
##'
##' KEGG2Sym('smu:SMU_20')
##'
##' KEGG2Sym('smu:SMU_408')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml2 read_xml xml_find_all xml_attr xml_text
##' @importFrom urltools url_encode
##' @importFrom magrittr %>%
##' @importFrom stringr str_detect str_trim
##' @importFrom KEGGAPI webTable
##' @export
##'
KEGG2Sym <- function(geneKEGG) {

  ## get id without KEGG species
  id <- geneKEGG %>%
    strsplit(split = ':') %>%
    sapply('[[', 2)

  ## get symbol
  KEGGsym <- geneKEGG %>%
    paste0('http://rest.kegg.jp/list/', .) %>%
    webTable(ncol = 2) %>%
    .[1, 2]

  ## split by ';' and ','
  ## if KEGGsym is '', then NULL will return. c('test', NULL) is equal to 'test'
  ele <- KEGGsym %>%
    strsplit(split = ';', fixed = TRUE) %>% ## split with ';'
    unlist %>%
    strsplit(split = ',', fixed = TRUE) %>% ## split with ','
    unlist %>%
    str_trim %>%
    .[!str_detect(., ' ')] %>% ## remove blank space
    c(id, .) ## add raw KEGGID

  return(ele)
}


##' KEGG symbol to Biocyc gene ID
##'
##' Try to convert the KEGG symbol to BioCyc gene ID.
##'
##' @title KEGG symbol to BioCyc ID.
##' @param symbol The gene symbol extract from KEGG.
##' @inheritParams ConvGeneKEGG2BioCyc
##' @return A \code{list}. The BioCyc gene ID or "0", if gene is not found.
##' @examples
##' Sym2BioCyc('dnaA', 'SMUT210007')
##'
##' Sym2BioCyc('SMU_06', 'SMUT210007')
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom xml2 read_xml xml_find_all xml_attr xml_text
##' @importFrom urltools url_encode
##' @importFrom magrittr %>%
##' @export
##'
Sym2BioCyc <- function(symbol, speBioCyc) {

  url <- paste0('http://websvc.biocyc.org/xmlquery?query=[x:x<-', speBioCyc, '^^genes,x^name%3D"', symbol, '"]&detail=full') %>%
    url_encode

  symbolxml <- read_xml(url)

  gene <- symbolxml %>%
    xml_find_all('//Gene[@ID][@frameid]') %>%
    xml_attr('frameid')

  protein <- symbolxml %>%
    xml_find_all('//Protein/@frameid') %>%
    xml_text

  res <- list(gene = gene,
              protein = protein,
              url = url)

  return(res)
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

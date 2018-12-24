##' Get node values.
##'
##' Get node values from the XML file.
##' @title Get all nodes value directly from the XML file.
##' @param xmlFile XML file.
##' @param nodePath The XPath of nodeset (one or mutiple nodes).
##' @return Nodeset value
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom XML getNodeSet xmlValue
##' @keywords internal
##' 
xmlNodeVal <- function(xmlFile, nodePath){

  nodeSet <- getNodeSet(xmlFile, nodePath)
  nodeValue <- sapply(nodeSet, xmlValue)

  return(nodeValue)

}

##' Get attribute values.
##'
##' Get attribute values from the XML file.
##' @title Get all nodes attributes directly from the XML file.
##' @param xmlFile XML file.
##' @param nodePath The XPath of nodeset (one or mutiple nodes).
##' @param attrName Attributes name
##' @return Nodeset attributes
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom XML getNodeSet xmlValue
##' @keywords internal
##' 
xmlNodeAttr <- function(xmlFile, nodePath, attrName){

  nodeSet <- getNodeSet(xmlFile, nodePath)
  nodeAttr <- sapply(nodeSet, xmlGetAttr, name = attrName)

  return(nodeAttr)

}


##' Test if the input vector's length is 0
##'
##' To test the length of input vector, if the length is 0, return "trueVal", else return "falseval"
##' @title Test length is 0 or not
##' @param inputVal vector
##' @param trueVal return this value, if the length of "inputVal" is 0.
##' @param falseVal return this value, if the length of "falseVal" is not 0.
##' @return "trueVal" or "falseVal"
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @keywords internal
##' 
testLen <- function(inputVal, trueVal, falseVal) {

  if (length(inputVal) == 0) {
    return(trueVal)
  } else {
    return(falseVal)
  }
}

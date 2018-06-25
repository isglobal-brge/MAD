#' Add annotation to a given object through the attribute 'gen.info'
#' 
#' @param x any object (usually 'parGADA' object)
#' @return the same object with an extra attribute called 'gen.info' including the annotation

addAnnot <- function(x){
  load(file.path(x, "SBL/gen.info.Rdata"))
  attr(x, "gen.info") <- gen.info
  return(x)
}

#' @title hlpFlattenDF
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @return SPATA object
#' @examples 
#' @export

hlpFlattenDF <- function(data.frame){
  map(.x=names(data.frame), .f=function(tab){data.frame[,tab]}) %>% unlist()
}


#' @title hlpList2Array
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @return SPATA object
#' @examples 
#' @export
#' 
hlpList2Array = function(input.list){  #input a list of lists
  rows.cols <- dim(input.list[[1]])
  sheets <- length(input.list)
  output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
  colnames(output.ary) <- colnames(input.list[[1]])
  row.names(output.ary) <- row.names(input.list[[1]])
  return(output.ary)    # output as a 3-D array
}
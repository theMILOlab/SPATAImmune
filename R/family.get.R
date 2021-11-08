#' @title GetVDJ
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @return SPATA object
#' @examples 
#' @export

GetVDJ <- function(object, constant=NULL){
  sample <- object@samples
  if(is.null(constant)){
    if(any(c(object@data[[sample]] %>% names())=="VDJ")){
      return(object@data[[sample]]$VDJ$IgBlast)
    }else{SPATAImmune::verbose.text("Data are not stored in object")}
  }else{
    object@data[[sample]]$VDJ$Clonality[[constant]]$IgBlast
  }


  
}

#' @title GetVDJ
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @return SPATA object
#' @examples 
#' @export
GetColors <-  function(n) {hues = seq(15, 375, length = n + 1);hcl(h = hues, l = 65, c = 100)[1:n]}

#' @title GetClonalityAnalysis
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @return SPATA object
#' @examples 
#' @export
GetClonalityAnalysis <- function(object, constant="TRB"){
  sample <- object@samples
  return(object@data[[sample]]$VDJ$ClusterValidation[[constant]] )
  
}


#' @title verbose.text
#' @author Dieter Henrik Heiland
#' @description Add on SPATAImmune
#' @inherit verbose prams 
#' @return 
#' @examples 
#' 
#' @export
#' 
#' 
#' 
verbose.text <- function(message){ message(paste0(Sys.time()," -- SPATA-Immune --  ", message)) }


#' @title check.SPATA
#' @author Dieter Henrik Heiland
#' @description Add on SPATAImmune
#' @inherit verbose prams 
#' @return 
#' @examples 
#' 
#' @export
#' 
#' 
#' 
check.SPATA <- function(){
  if("SPATA2" %in% rownames(installed.packages())){
    verbose.text("SPATA2 was found")
  }else{ stop(verbose.text("SPATA2 is required for SPATAImmune, please install teh SPATA2 package: https://github.com/theMILOlab/SPATA2 "))
    }
}



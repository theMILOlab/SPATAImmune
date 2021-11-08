#' @title readsamfiles
#' @author Dieter Henrik Heiland
#' @description Add on SPATAImmune
#' @inherit verbose prams 
#' @return 
#' @examples 
#' 
#' @export
#' 
readsamfiles<- function(sam.file){
  SPATAImmune::check.SPATA()
  
  SPATAImmune::verbose.text("The import function was adopted from the https://github.com/VladaMilch/pdProbeRemap package, please cite")
  
  line2vals<- function(x,nfields=11){
    x = base::strsplit(x,"\t",fixed=TRUE)[[1]]
    c(x[1:nfields],base::paste(x[-c(1:nfields)],collapse="\t"))
  }
  
  
  SPATAImmune::verbose.text("reading file...")
  x = base::readLines(sam.file)
  SPATAImmune::verbose.text("fetching header...")
  headerpos = base::grep("^@",x)
  header = x[headerpos]
  SPATAImmune::verbose.text("converting header...")
  
  header = list("HD" = base::lapply(base::gsub("^@HD\t","",header[base::grep("^@HD",header)]),function(x)base::strsplit(x,"\t")[[1]]),
                "SQ" = base::lapply(base::gsub("^@SQ\t","",header[base::grep("^@SQ",header)]),function(x)base::strsplit(x,"\t")[[1]]),
                "RG" = base::lapply(base::gsub("^@RG\t","",header[base::grep("^@RG",header)]),function(x)base::strsplit(x,"\t")[[1]]),
                "PG" = base::lapply(base::gsub("^@PG\t","",header[base::grep("^@PG",header)]),function(x)base::strsplit(x,"\t")[[1]]),
                "CO" = base::gsub("^@CO\t","",header[base::grep("^@CO",header)])
  )
  SPATAImmune::verbose.text("fetching data...")
  x = base::lapply(x[-headerpos],line2vals)
  
  SPATAImmune::verbose.text("convert to data.frame...")
  x = matrix(base::unlist(x),nrow=length(x),ncol=length(x[[1]]),byrow=TRUE,dimnames=list(NULL,c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","TLEN","SEQ","QUAL","ATTR")))
  x = as.data.frame(x,stringsAsFactors=FALSE)
  x$FLAG<- base::as.integer(x$FLAG)
  x$POS<- base::as.integer(x$POS)
  x$MAPQ<- base::as.integer(x$MAPQ)
  x$PNEXT<- base::as.integer(x$PNEXT)
  x$TLEN<- base::as.integer(x$TLEN)
  SPATAImmune::verbose.text("done!\n")
  base::return(list("header"=header,"x" = x))
}


#' @title create.IAF
#' @author Dieter Henrik Heiland
#' @description This function creates a TCR aligned file from the alignment (spatial bc) and the TCR-sequences
#' @inherit verbose prams 
#' @return 
#' @examples 
#' 
#' @export
#'
#'
create.IAF <- function(Binner.file, sam.file){
  SPATAImmune::check.SPATA()
  SPATAImmune::verbose.text("Check sam files")
  
  #Check sam files
  if(sam.file %>% class() == "list"){SPATAImmune::verbose.text("sam file was loaded") }else{
    if(!base::file.exists(sam.file)) stop(SPATAImmune::verbose.text("No valid sam file found"))
    SPATAImmune::readsamfiles(sam.file)
  }
  
  #Check Binner.file
  if(is.null(Binner.file)) stop(SPATAImmune::verbose.text("Binner file not determined"))
  if(!base::file.exists(Binner.file)) stop(SPATAImmune::verbose.text("Binner file not valid"))
  #Read file
  Binner.file <- read.txt(Binner.file)
  
  
  
  
  
}



#' @title import.sam.files
#' @author Dieter Henrik Heiland
#' @description Add on SPATAImmune
#' @inherit verbose prams 
#' @return 
#' @examples 
#' 
#' @export
#'
#'
import.sam.files <- function(object, IAF.file){}




#' @title addReArrangmentMatrix
#' @author Dieter Henrik Heiland
#' @description This function transformes a input to a expression matrix and import it to the spata object
#' @param vdj.df The data frame with barcodes and the ReArrangment parameters
#' @param barcodes Colname of the barcodes 
#' @param target.col Colname of the target enrichment colum
#' @param mtr.name The name of the importet dataset
#' @inherit verbose prams 
#' @return SPATA object
#' @examples 
#' 
#' @export
#'
#'
addReArrangmentMatrix <- function(object, 
                                  vdj.df,
                                  target.col="reAr",
                                  barcodes="BC",
                                  mtr.name=NULL,
                                  set.active=T){
  
  #Check input
  if(!ncol(vdj.df[,c(barcodes, target.col )])==2) stop(" Input vdj colnames corrupt")
  if(is.null(mtr.name)){mtr.name="No_Defined"}
  
  #Reshape
  acast.df=vdj.df %>% dplyr::select({{barcodes}}, {{target.col}})
  names(acast.df) <- c("barcodes", "target.col")
  mat.out <- reshape2::acast(barcodes ~ target.col, data= acast.df) %>% t()
  
  # Add missing rows (BCs)
  bc.all <- SPATA2::getCoordsDf(object) %>% dplyr::pull(barcodes)
  BC.missing <- bc.all[! bc.all %in%  colnames(mat.out)] 
  missing <- matrix(0, nrow(mat.out), length(BC.missing) )
  colnames(missing) <- BC.missing
  mat.out <- cbind(mat.out, missing)
  object <- SPATA2::addExpressionMatrix(object, expr_mtr=mat.out, mtr_name=mtr.name)
  if(set.active==T){
    object <- setActiveExpressionMatrix(object, mtr_name=mtr.name)
  }
  
}



#' @title ImportVDJ -Deprecated - 
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @param object The data frame with barcodes and the ReArrangment parameters
#' @param barcodes.df Data frame of the barcodes, with the following cols: c("ID", "barcodes"): ID: Sequencing ID, barcodes: The annotated Barcode
#' @param IgBlast.df Optimized IgBlast output with the following cols: c("ID", "c", "v", "d", "j", "CDR3"). ID: Sequencing ID, c: Constant Region, v,d,j calls, and the CDR3 sequence
#' @param Pychoper If the Pychopper option was used before IgBlast (in the pipline)
#' @inherit verbose prams 
#' @return SPATA object
#' @examples 
#' 
#' @export



ImportVDJ <- function(object,
                      barcodes.df,
                      IgBlast.df,
                      Pychoper=F,
                      Filter.cdr3=T,
                      ram=8,
                      cors=16){
  
  SPATA2::check_object(object)
  
  #Check barcodes.df
  if(!ncol(barcodes.df[,c("ID", "barcodes" )])==2) stop(" Input barcodes.df colnames corrupt")
  if(!ncol(IgBlast.df[,c("ID", "c", "v", "d", "j", "CDR3")])==6) stop(" Input IgBlast.df colnames corrupt")
  
  SPATAImmune::verbose.text("Input ok... ")
  
  
  SPATAImmune::verbose.text("Align TCR and Barcodes ")
  vdj <- IgBlast.df
  
  if(Filter.cdr3==T){
    SPATAImmune::verbose.text("Remove non CDR3_aa Reads ")
    vdj <- 
      vdj %>% 
      dplyr::filter(CDR3!="")
  }
  
  
  #After Pychop the header needs to be optimized
  if(Pychoper==T){
    SPATAImmune::verbose.text(" Manipulate header if Pychopper was used ")
    
    #Clean Seq ID from Pychoper manipulation
    ID <- 
      vdj %>% 
      dplyr::select(ID) %>% 
      tidyr::separate(., col="ID", into=c("ID", "e"), sep="runid") %>% 
      tidyr::separate(., col="ID", into=c("f", "ID"), sep="[|]") %>% 
      dplyr::pull(ID)
    
    
    #options(future.fork.enable = TRUE)
    #future::plan("multiprocess", workers = cors)
    #options(future.globals.maxSize = ram * 1024^2)
    #message("... Run multicore ... ")
    
    #ID <- furrr::future_map(.x=ID, .f=function(x){
    #  out <- stringr::str_split(x, pattern = "runid=") %>% unlist()
    #  out <- stringr::str_split(out[1], pattern = "[|]") %>% unlist()
    #  return(out[2])
    #}, .progress=T)

    vdj$ID <- ID
  }
  
  
  SPATAImmune::verbose.text("Load Barcodes from SPATA object ")
  all.bc <- object %>% SPATA2::getCoordsDf()
  
  vdj <- 
    vdj %>% 
    dplyr::left_join(., barcodes.df, by="ID") %>% 
    dplyr::mutate(barcodes=paste0(barcodes, "-1")) %>% 
    dplyr::filter(barcodes %in% all.bc$barcodes)

  SPATAImmune::verbose.text("Clean data and remove non-full Arrangments ")
  
  vdj <- 
    vdj %>% 
    dplyr::mutate(Arrangment=dplyr::case_when(
      c=="TRA" & v!="" & j!=""  ~ "1",
      c=="TRB" & v!="" & d!=""& j!=""  ~ "1",
      c=="TRD" & v!="" & j!=""  ~ "1",
      c=="TRG" & v!="" & d!=""& j!=""  ~ "1",
      TRUE ~ "0"
    )) %>% 
    filter(Arrangment==1)
  
  
  SPATAImmune::verbose.text("Add AA-Seq ")
  
  vdj$CDR3_aa <-  
    Biostrings::translate(vdj$CDR3 %>% Biostrings::DNAStringSet()) %>% 
    as.data.frame() %>% 
    dplyr::pull(x)
  
  
  SPATAImmune::verbose.text("Add data to object")
  
  sample <- object@samples
  
  object@data[[sample]]$VDJ <- list()
  object@data[[sample]]$VDJ$IgBlast <- vdj
  
  return(object)
  
  
}



#' @title ImportSPTCRseq
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @param object The data frame with barcodes and the ReArrangment parameters
#' @param path The path to the SPTCR-seq pipeline output


#' @return SPATA object
#' @examples 
#' 
#' @export


ImportSPTCRseq <- function(object,
                           sample_name,
                           path){
  
  SPATA2::check_object(object)
  
  SPATAImmune::verbose.text("Load SPTCR seq Pipeline data")
  vdj <- read.csv(paste0(path,"/ClusterCorrect/",sample_name,"_CORRECTED_umi_corrected_count_table.csv"))
  vdj <- vdj %>% mutate(barcodes=paste0(vdj$Spatial.Barcode, "-1"))
  
  
  SPATAImmune::verbose.text("Check Barcodes")
  
  bc <- object %>% SPATA2:::getBarcodes()
  barcode.stat <- table(vdj$barcodes %in% bc) %>% as.data.frame()
  message <- paste0(round(c(barcode.stat$Freq[2]/barcode.stat$Freq[1]*100),digits = 2) , " % of the barcodes found in the stRNA-seq object")
  SPATAImmune::verbose.text(message)
  
  vdj_filtered <- 
    vdj %>% 
    dplyr::filter(barcodes %in% bc) %>% 
    dplyr::select(barcodes, Locus,V,D,J,CDR3_aa, UMI.Corrected)
  
  
  #vdj_filtered %>% group_by(barcodes) %>% summarise(sum=length(barcodes)) %>% arrange(desc(sum))
  #vdj_filtered %>% group_by(CDR3_aa) %>% summarise(sum=length(CDR3_aa)) %>% arrange(desc(sum))
  
  names(vdj_filtered) <- c("barcodes", "c", "v", "d", "j", "CDR3_aa", "Count")
  
  SPATAImmune::verbose.text("Add data to object")
  sample <- object@samples
  object@data[[sample]]$VDJ <- list()
  object@data[[sample]]$VDJ$SPTCR <- vdj_filtered
  
  return(object)
  
  
}




































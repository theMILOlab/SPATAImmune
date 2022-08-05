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

#' @title GetSPTCR
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline 
#' @return SPATA object
#' @examples 
#' @export

GetSPTCR <- function(object){
  sample <- object@samples
    if(any(c(object@data[[sample]] %>% names())=="VDJ")){
      return(object@data[[sample]]$VDJ$SPTCR)}else{SPATAImmune::verbose.text("Data are not stored in object")}
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
GetClonalityAnalysis <- function(object, 
                                 constant="TRB"){
  sample <- object@samples
  return(object@data[[sample]]$VDJ$ClusterValidation[[constant]] )
  
}





#' @title getSpatialRegression
#' @author Dieter Henrik Heiland
#' @param object SPATA object.
#' @param features features to perform spatial regression.
#' @param model The model for spatial weighted regression. classical:Pearson (Not recommended), CCA model, Spatial Lag model: lmSLX, dist 
#' @param smooth Kernel.
#' @param normalize Normalize Input (Not recommended).
#' @description Wrapper function for spatial weighted regression analysis
#' @return Correlation matrix of NxN input
#' @export
#' 
getSpatialRegression <- function(object, 
                                 features, 
                                 model, 
                                 smooth=F, 
                                 normalize=F){
  
  #1. transform data into sp format
  color_var <- SPATA2::hlpr_join_with_aes(object, df= SPATA2::getCoordsDf(object), color_by = features, normalize=normalize, smooth=smooth)
  coords <- color_var[,c("x", "y")] %>% as.data.frame()
  data <- color_var[,features] %>% as.data.frame()
  rownames(coords)=rownames(data) <- color_var$barcodes
  sp_obj <- sp::SpatialPointsDataFrame(coords = coords, data=data)
  
  
  #2. Get Neighbors
  nc <- spdep::knearneigh(sp_obj, k=6, longlat = T)
  nc.nb <- spdep::knn2nb(nc)
  nb.list <- spdep::nb2listw(nc.nb)
  nb.mat <- spdep::nb2mat(nc.nb)
  runCCADist <- function(object,color_var, a,b){
    sample <- SPATA2::getSampleNames(object)
    data <- color_var[,c(a,b)] %>% as.data.frame()
    names(data) <- c("a", "b")
    rownames(data) <- color_var$barcodes
    
    coord_spots <- object %>% SPATA2::getCoordsDf()
    coord_spots <- coord_spots[, c("row", "col")] %>% as.data.frame()
    rownames(coord_spots) <- object %>% SPATA2::getCoordsDf() %>% pull(barcodes)
    
    coord_spots <- cbind(coord_spots, data[rownames(coord_spots), ])
    coord_spots2 <- coord_spots
    coord_spots2$a <- coord_spots2$b*c(-1)
    
    getKernel <- function(coord_spots, var){
      mat <- reshape2::acast(row~col, data=coord_spots, value.var = var)
      mat[is.na(mat)] <- 0
      mat %>%  scales::rescale(., c(0,1)) %>% oce::matrixSmooth() 
    }
    
    
    
    real <- proxy::dist(x=getKernel(coord_spots, "a"),
                        y=getKernel(coord_spots, "b"))
    
    X <- getKernel(coord_spots, "a")
    Y <- getKernel(coord_spots, "b")
    cca <- CCA::matcor(X, Y)
    
    return(c(cca$XYcor[!is.na(cca$XYcor)] %>% mean(), mean(na.omit(as.vector(real)))))
    
    
  }
  
  #classical regression model
  if(model=="classical"){
    
    mat.cor <- cor(data)
    
    return(mat.cor)
    
  }
  
  if(model=="lmSLX"){
    
    mat.cor <- matrix(NA, length(features), length(features))
    rownames(mat.cor) = colnames(mat.cor) <- features
    
    for(i in 1:length(features)){
      for(j in 1:length(features)){
        
        first_feat <- features[i]
        second_feat <- features[j]
        if(first_feat==second_feat){mat.cor[first_feat,second_feat]=0}else{
          formula <- as.formula(paste(first_feat,"~", second_feat))
          reg1 <- lm(formula, data = sp_obj)
          
          reg2 = spatialreg::lmSLX(reg1, data = sp_obj, nb.list)
          sum_mod_2 <- summary(reg2)
          cor_lag2 <- sum_mod_2$coefficients[3,1]
          mat.cor[first_feat,second_feat]=as.numeric(cor_lag2)}
      }
    }
    
    return(mat.cor)
    
  }
  
  if(model=="lagsarlm"){
    if(length(features)>=2) stop("Run time to long, reduce number of features")
    
    first_feat <- features[1]
    second_feat <- paste(features[2:length(features)], collapse = "+")
    formula <- as.formula(paste(first_feat,"~", second_feat))
    reg1 <- lm(formula, data = sp_obj)
    sum_mod_1 <- summary(reg1)
    cor_lag1 <- sum_mod_1$coefficients[2,]
    
    #lagsarlm
    reg3 = spatialreg::lagsarlm(reg1, data = sp_obj, nb.list)
    sum_mod_3 <- summary(reg3)
    mat.cor <- data.frame(estimate=as.numeric(sum_mod_3$coefficients[2]), 
                          p.value=as.numeric(sum_mod_3$Wald1$p.value))
    
  }
  
  if(model=="errorsarlm"){
    if(length(features)>=2) stop("Run time to long, reduce number of features")
    
    first_feat <- features[1]
    second_feat <- paste(features[2:length(features)], collapse = "+")
    formula <- as.formula(paste(first_feat,"~", second_feat))
    reg1 <- lm(formula, data = sp_obj)
    sum_mod_1 <- summary(reg1)
    cor_lag1 <- sum_mod_1$coefficients[2,]
    
    #errorsarlm
    reg4 = spatialreg::errorsarlm(reg1, data = sp_obj, nb.list)
    sum_mod_4 <- summary(reg4)
    mat.cor <- data.frame(estimate=as.numeric(sum_mod_4$coefficients[2]), 
                          p.value=as.numeric(sum_mod_4$Wald1$p.value))
    
  }
  
  if(model=="CCA"){
    
    #CCA model
    mat.cor <- matrix(NA, length(features), length(features))
    rownames(mat.cor) = colnames(mat.cor) <- features
    for(i in features){
      for(j in features){
        if(i==j){mat.cor[i,j] <- 1}else{mat.cor[i,j] <- runCCADist(object, a=i, b=j, color_var=color_var)[1]}
      }
    }
    
    return(mat.cor)
    
  }
  
  if(model=="dist"){
    
    mat.cor <- matrix(NA, length(features), length(features))
    rownames(mat.cor) = colnames(mat.cor) <- features
    for(i in features){
      for(j in features){
        if(i==j){mat.cor[i,j] <- 0}else{mat.cor[i,j] <- runCCADist(object, a=i, b=j, color_var=color_var)[2]}
      }
    }
    
    return(mat.cor)
    
  }
  
  return(mat.cor)
  
  
}









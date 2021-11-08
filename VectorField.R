library(metR)
library(ggplot2)
library(data.table)
library(tidyverse)
library(SPATA2)


# Real spatial dataset
object <- readRDS("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions/275_T_SPATA_CNV_Pred.RDS")


# get VF
parameter= "HM_G2M_CHECKPOINT"
#parameter = "HM_HYPOXIA"
df <- joinWithGeneSets(object, gene_sets = parameter)

plot.VF.spatial <- function(df, 
                            object,
                            parameter, 
                            dist.spot=10,
                            color.by,
                            vector.plot=F, 
                            stream.plot=T){

#Define color parameter
color.points <- df %>% pull(!!sym(color.by))

df <- df %>% as.data.frame() %>% dplyr::select(barcodes,x,y,{{parameter}})


get.vector.field <- function(object, df, parameter,dist.spot=dist.spot){
  
  

# functions ---------------------------------------------------------------

  # Get NN
  getSurroundedSpots <- function(object){
    of_sample <- SPATA2::getSampleNames(object)
    coords <- SPATA2::getCoordsDf(object)
    bc_origin <- coords$barcodes
    bc_destination <- coords$barcodes
    
    # get grouped data.frame with all barcode combinations
    
    dfgr <-
      tidyr::expand_grid(bc_origin, bc_destination) %>%
      dplyr::left_join(x = ., y = dplyr::select(coords, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>%
      dplyr::left_join(x = ., y = dplyr::select(coords, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>%
      dplyr::mutate(distance = base::round(base::sqrt((xd - xo)^2 + (yd - yo)^2), digits = 0)) %>%
      dplyr::group_by(bc_origin)
    
    # filter barcodes that are surrounded by a total of six spots
    
    sufficient_bc <-
      dplyr::slice_min(.data = dfgr, order_by = distance, n = 7) %>%
      dplyr::group_by(bc_origin, distance) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::filter(distance != 0 & count == 6) %>%
      dplyr::pull(bc_origin)
    
    # filter for barcodes 
    
    final_df <-
      dplyr::filter(dfgr, bc_origin %in% sufficient_bc) %>%
      dplyr::slice_min(order_by = distance, n = 7) %>%
      dplyr::mutate(sample = {{of_sample}}) %>%
      dplyr::select(sample, dplyr::everything()) 
    
    base::return(final_df)
    
  }
  

# Prepare data ------------------------------------------------------------
  NN.file <- getSurroundedSpots(object)
  library(furrr)
  options(future.fork.enable = TRUE)
  plan("multiprocess", workers = 8)
  supportsMulticore()
  options(future.globals.maxSize = 50000 * 1024^2)
  message("... Run multicore ... ")
  
  VF <- furrr::future_map(.x=1:nrow(df), .f=function(i){
    
    bc <- df[i, c("barcodes")]
    cc <- df[i, c("x", "y")]
    
    #NN <- NN.file %>% filter(bc_origin=={{bc}}) %>% pull(bc_destination)
    
    NN <- NN.file %>% 
      filter(xo < cc$x+dist.spot & xo > cc$x-dist.spot) %>% 
      filter(yo < cc$y+dist.spot & yo > cc$y-dist.spot) %>%
      pull(bc_destination)
    
    
    NN.df <- df %>% filter(barcodes %in% NN) %>% as.data.frame()
    
    # create vector
    V <- -c(
      as.numeric(cc) - c(NN.df$x[which.max(NN.df[,parameter])], NN.df$y[which.max(NN.df[,parameter])])
    )
    
    #V <- c(NN$x[which.max(NN$z)], NN$y[which.max(NN$z)]) * cor[i, c("z")]
    if(length(V)==0){
      out <- data.frame(barcodes=bc, t.x=0, t.y=0)
    }else{out <- data.frame(barcodes=bc, t.x=V[1], t.y=V[2])}
    
    
    return(out)
    
    
    
  }, .progress = T) %>% 
    do.call(rbind, .) %>% 
    as.data.frame()
  
  
  out <- cbind(df, VF)
  out[is.na(out)] <- 0
  
  return(out)
  
  
  
  
}
VF <- get.vector.field(object, df, parameter=parameter, dist.spot=dist.spot)

message( "
         ... Create Plot ...
         ")
VF <- VF %>% dplyr::select(x,y,{{parameter}}, t.x, t.y) %>% rename("parameter":=!!sym(parameter))
if(vector.plot==T){
  ggplot(data=VF, aes(x,y))+
    geom_point(data=VF , mapping=aes(x,y, color=color.points), size=3, alpha=0.5)+
    theme_classic()+
    geom_vector(aes(dx = t.x, dy = t.y)) +
    scale_mag()
}
if(stream.plot==T){

plot.streamlines <- function(VF,size.arrow=1){
  
  drifter.split.sf = VF %>% 
    sf::st_as_sf(coords = c("x", "y"))
  
  
  drifter.grid = drifter.split.sf %>% 
    st_make_grid(n = c(70,60))%>%
    st_sf()
  
  drifter.split.sf.se = drifter.split.sf%>% filter(parameter!=0)
  
  drifter.gridded = drifter.grid %>% 
    mutate(id = 1:n(), contained = lapply(st_contains(st_sf(geometry),drifter.split.sf.se),identity),
           obs = sapply(contained, length),
           u = sapply(contained, function(x) {median(drifter.split.sf.se[x,]$t.x, na.rm = TRUE)}),
           v = sapply(contained, function(x) {median(drifter.split.sf.se[x,]$t.y, na.rm = TRUE)})) 
  
  
  
  drifter.gridded = drifter.gridded %>% dplyr::select(obs, u, v) %>% na.omit()
  
  ## obtain the centroid coordinates from the grid as table
  coordinates = drifter.gridded %>% 
    st_centroid() %>% 
    st_coordinates() %>% 
    as_tibble() %>% 
    rename(x = X, y = Y)
  
  ## remove the geometry from the simple feature of gridded drifter dataset
  st_geometry(drifter.gridded) = NULL
  
  ## stitch together the extracted coordinates and drifter information int a single table for SE monsoon season
  current.gridded.se = coordinates %>% 
    bind_cols(drifter.gridded) %>% 
    mutate(season = "SE")
  
  ## bind the gridded table for SE and NE
  ## Note that similar NE follow similar procedure, hence not shown in the post
  drifter.current.gridded = current.gridded.se %>% 
    bind_rows(current.gridded.se)
  
  
  
  
  ## select grids for SE season only
  drf.se = drifter.current.gridded %>% filter(season == "SE")
  
  ## interpolate the U component
  u.se = oce::interpBarnes(x = drf.se$x, y = drf.se$y, z = drf.se$u, gamma=0.5)
  
  ## obtain dimension that determine the width (ncol) and length (nrow) for tranforming wide into long format table
  dimension = data.frame(lon = u.se$xg, u.se$zg) %>% dim()
  
  ## make a U component data table from interpolated matrix
  u.tb = data.frame(lon = u.se$xg, 
                    u.se$zg) %>% 
    gather(key = "lata", value = "u", 2:dimension[2]) %>% 
    mutate(lat = rep(u.se$yg, each = dimension[1])) %>% 
    dplyr::select(lon,lat, u) %>% as.tibble()
  
  ## interpolate the V component
  v.se = oce::interpBarnes(x = drf.se$x, 
                           y = drf.se$y, 
                           z = drf.se$v,
                           gamma=0.5)
  
  ## make the V component data table from interpolated matrix
  v.tb = data.frame(lon = v.se$xg, v.se$zg) %>% 
    gather(key = "lata", value = "v", 2:dimension[2]) %>% 
    mutate(lat = rep(v.se$yg, each = dimension[1])) %>% 
    dplyr::select(lon,lat, v) %>% 
    as.tibble()
  
  ## stitch now the V component intot the U data table and compute the velocity
  uv.se = u.tb %>% 
    bind_cols(v.tb %>% dplyr::select(v)) %>% 
    mutate(vel = sqrt(u^2+v^2))
  
  library(oce)
  
  if(color.points %>% class()=="factor"){
    ggplot() +
      geom_point(data=VF, mapping=aes(x,y, color=color.points), size=5, alpha=0.5)+
      metR::geom_streamline(data = uv.se, 
                            aes(x = lon, y = lat, dx = u, dy = v),
                            size=size.arrow,
                            arrow.length = 0.5,
                            arrow.angle = 25,
                            arrow.type = "closed",
                            L = 3, res =0.5,
                            lineend = "round")+
      theme_void()
  }else{
    ggplot() +
      geom_point(data=VF, mapping=aes(x,y, color=color.points), size=5, alpha=0.5)+
      scale_color_viridis_c(guide = "none") +
      metR::geom_streamline(data = uv.se, 
                            aes(x = lon, y = lat, dx = u, dy = v),
                            size=size.arrow,
                            arrow.length = 0.5,
                            arrow.angle = 25,
                            arrow.type = "closed",
                            L = 3, res =0.5,
                            lineend = "round")+
      theme_void()
  }
  
  
  
  
}

plot.streamlines(VF,size.arrow=1)
}

}






#Not used:
get.spatial.grid <- function(path){
  
  grid.info <- paste0(
    path, "/outs/spatial/tissue_positions_list.csv"
  )
  
  if(file.exists(grid.info)){
    tl <- read.csv(grid.info, header=F)
  } else{ stop("no file")}
  
  scale.info <- paste0(
    path, "/outs/spatial/scalefactors_json.json"
  )
  
  if(file.exists(scale.info)){
  s.i <- suppressWarnings(rjson::fromJSON(paste(readLines(scale.info), collapse="")))
  } else{ stop("no file")}
  
  scale.f <- s.i$tissue_lowres_scalef
  
  ## scale data
  tl$V5 <- tl$V5*scale.f
  tl$V6 <- tl$V6*scale.f
  
  return(tl)
}



parameter= "HM_HYPOXIA"
df <- joinWithGeneSets(object, gene_sets = parameter)

plot.VF.spatial(df, parameter, spatial.grid, dist.spot=100)

parameter= "BP.GO_CEREBRAL_CORTEX_RADIALLY_ORIENTED_CELL_MIGRATION"
SPATA2::getGeneSets(object, index = "CELL_MIGRATION")
df <- joinWithGeneSets(object, gene_sets = parameter) %>% 
  left_join(., 
            joinWithFeatures(object, features="seurat_clusters") %>% dplyr::select(barcodes, seurat_clusters),
            by="barcodes")


plot.VF.spatial(df,object, parameter, dist.spot=100, color.by="seurat_clusters")




# Real spatial dataset
object <- readRDS("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions/313_T_SPATA_CNV_Pred.RDS")



parameter= "HM_HYPOXIA"
SPATA2::getGeneSets(object, index = "CELL_MIGRATION")
df <- joinWithGeneSets(object, gene_sets = parameter) %>% 
  left_join(., 
            joinWithFeatures(object, features="seurat_clusters") %>% dplyr::select(barcodes, seurat_clusters),
            by="barcodes")


plot.VF.spatial(df,object, parameter, dist.spot=100, color.by="seurat_clusters")

parameter= "BP.GO_CEREBRAL_CORTEX_RADIALLY_ORIENTED_CELL_MIGRATION"
SPATA2::getGeneSets(object, index = "CELL_MIGRATION")
df <- joinWithGeneSets(object, gene_sets = parameter) %>% 
  left_join(., 
            joinWithFeatures(object, features="seurat_clusters") %>% dplyr::select(barcodes, seurat_clusters),
            by="barcodes")


plot.VF.spatial(df , object, parameter, dist.spot=100, color.by="seurat_clusters")


object %>% SPATA2::plotSurfaceInteractive()



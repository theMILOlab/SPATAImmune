
# Plot functions ----------------------------------------------------------




# basic statistics --------------------------------------------------------

#' @title plot.donut
#' @author Dieter Henrik Heiland
#' @description SPATA-Immuno
#' @return 
#' @examples 
#' 
#' @export
#' 
plotDonut <- function(table){
  
  check.SPATA()
  verbose.text("Create Donut Plot")
  
  # Prepare data 
  data <- data.frame(table)
  base::names(data) <- c("category","count")
  data$fraction = data$count / base::sum(data$count)
  data$ymax = base::cumsum(data$fraction)
  data$ymin = c(0, utils::head(data$ymax, n=-1))
  
  # Make the plot
  ggplot2::ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    ggplot2::geom_rect() +
    ggplot2::coord_polar(theta="y") +
    ggplot2::xlim(c(2, 4))+ 
    ggplot2::theme_void()
  
  
  
}



#' @title River Plot
#' @author Dieter Henrik Heiland
#' @description SPATA-Immuno
#' @return 
#' @examples 
#' 
#' @export
#' 
plotRiver <- function(data, var1, var2, color1=NULL, color2=NULL) {
  
  require(dplyr)          # Needed for the count function
  require(riverplot)      # Does all the real work
  require(RColorBrewer)   # To assign nice colours
  
  
  
  names1 <- paste0("A_",levels(data %>% pull(!!sym(var1)) %>% as.character() %>% as.factor()))
  names2 <- paste0("B_",levels(data %>% pull(!!sym(var2)) %>% as.character() %>% as.factor()))
  
  #names1 %>% length()
  #names2 %>% length()
  
  #Check colors
  if(!is.null(color1)){if(length(color1)!=length(names1)) stop("color1 and levels do not match")}
  if(!is.null(color2)){if(length(color2)!=length(names2)) stop("color2 and levels do not match")}
  
  
  var1   <- as.numeric(data[, var1] %>% as.factor())
  #var1 %>% unique() %>% length()
  var2   <- as.numeric(data[, var2] %>% as.factor())
  #var2 %>% unique() %>% length()
  
  
  
  edges  <- data.frame(as.factor(var1), as.factor(var2 + max(var1, na.rm = T)))
  colnames(edges) <- c("N1", "N2")
  edges  <- edges %>% group_by(N1,N2) %>% count()
  colnames(edges) <- c("N1", "N2", "Value")
  
  
  nodes <- data.frame(
    ID     = c(1:(max(var1, na.rm = T) + max(var2, na.rm = T))),  
    x      =  c(rep(1, times = max(var1, na.rm = T)), 
                rep(2, times = max(var2, na.rm = T))),       
    labels = c(names1, names2) , 
    col    = c(if(!is.null(color1)){color1}else{colorRampPalette(brewer.pal(9, "Greys"))(max(var1, na.rm = T))} , 
               if(!is.null(color2)){color2}else{colorRampPalette(brewer.pal(9, "Reds"))(max(var2, na.rm = T))} ) ,
    stringsAsFactors = FALSE)
  
  if(is.null(color1)){nodes$col <- paste(nodes$col, 95, sep = "")}
  
  ds <- default.style()
  ds$srt="0"
  river <- makeRiver(nodes, edges %>% as.data.frame())
  
  return(plot(river, default_style=ds))
}

#' @title plotQualityIgBlast
#' @author Dieter Henrik Heiland
#' @description SPATA-Immuno
#' @return 
#' @examples 
#' 
#' @export
#' 
plotQualityIgBlast <- function(vdj){
  
  require(patchwork)
  
  #Raw non spatial alligned locus reads
  p1 <- ggplot(data=vdj %>% filter(c!=""), aes(x=c, fill=c)) + geom_bar() +theme_classic()+ggtitle("Numver of Reads per Constant Region")
  
  #Raw non spatial alligned Full CDR3 seq
  res_raw <- vdj %>% filter(c!="") %>%  mutate(Detected=ifelse(CDR3!="", "1", "0")) %>% group_by(c) %>% count(Detected)
  p2 <- ggplot(data=res_raw, aes(x=c, y=n, fill=Detected)) + geom_bar(position = "fill", stat = "identity") + theme_classic()+ggtitle("Full CDR3 Region sequenced")
  
  #Raw non spatial alligned Full VDJ Arrangment
  res_raw <- 
    vdj %>% 
    filter(c!="") %>% 
    dplyr::mutate(number=1:nrow(.)) %>% 
    dplyr::mutate(Arrangment=dplyr::case_when(
      c=="TRA" & v!="" & j!=""  ~ "1",
      c=="TRB" & v!="" & d!=""& j!=""  ~ "1",
      c=="TRD" & v!="" & j!=""  ~ "1",
      c=="TRG" & v!="" & d!=""& j!=""  ~ "1",
      TRUE ~ "0"
    )) %>% 
    group_by(c) %>% 
    count(Arrangment)
  p3=ggplot(data=res_raw, aes(x=c, y=n, fill=Arrangment)) + geom_bar(position = "fill", stat = "identity") + theme_classic()+ggtitle("Full TCR Arrangment")
  
  
  #Quality 
  p <- p1+p2+p3
  
  return(p)
  
  
}



#' @title plotVDJCirclize
#' @author Dieter Henrik Heiland
#' @description SPATA-Immuno
#' @return 
#' @examples 
#' 
#' @export
#' 
plotVDJCirclize <- function(object, 
                            select.c="TRB", 
                            arrange.by=NULL, 
                            top=200, 
                            output=getwd()){
  
  data <- SPATAImmune::GetVDJ(object)
  if(is.null(data)){data <- SPATAImmune::GetSPTCR(object)}
  if(is.null(data)) stop ("No VDJ or SPTCR slot in object")
  
  # If required select spots
  if(!is.null(arrange.by)){
    bc.join <- joinWithFeatures(object, features = arrange.by) %>% select(barcodes, {{arrange.by}})
    data <- 
      data %>% 
      dplyr::left_join(.,bc.join, by="barcodes") %>% 
      as.data.frame()
    
    purrr::map(.x=unique(data[,arrange.by]),
               .f=function(select){
                 data <- data %>% dplyr::filter(!!sym(arrange.by)=={{select}})
                 
                 if(select.c == "TRB" | select.c == "TRD"){
                   data <- 
                     data %>% 
                     dplyr::filter(c==select.c) %>% 
                     dplyr::select(c,v,d,j, barcodes) %>% 
                     dplyr::mutate(reAr=paste0(v,"_",d,"_",j))
                   
                   if(nrow(data)>top){
                     
                     select.top <- 
                       data %>% 
                       dplyr::group_by(reAr) %>% 
                       dplyr::count() %>% 
                       dplyr::arrange(desc(n)) %>% 
                       head(top) %>% 
                       dplyr::pull(reAr)
                     
                     data <- 
                       data %>% 
                       dplyr::filter(reAr %in% select.top)}
                   
                 }else{
                   
                   data <- 
                     data %>% 
                     dplyr::filter(c==select.c) %>% 
                     dplyr::select(c,v,j,barcodes) %>% 
                     dplyr::mutate(reAr=paste0(v,"_",j))
                   
                   if(nrow(data)>top){
                     
                     select.top <- 
                       data %>% 
                       dplyr::group_by(reAr) %>% 
                       dplyr::count() %>% 
                       dplyr::arrange(desc(n)) %>% 
                       head(top) %>% 
                       dplyr::pull(reAr)
                     
                     data <- 
                       data %>% 
                       dplyr::filter(reAr %in% select.top)
                     
                     
                   }
                   
                 }
                 
                 
                 
                 #save plot
                 setwd(output)
                 pdf(paste0("VDJ_Plot_",arrange.by,"_",select,".pdf"), useDingbats=F)
                 
                 order_v <- data %>% group_by(v) %>% summarise(sum=length(v)) %>% arrange(desc(sum)) %>% ungroup() %>% pull(v)
                 order_j <- data %>% group_by(j) %>% summarise(sum=length(j)) %>% arrange(desc(sum)) %>% ungroup() %>% pull(j)
                 
                 mat <- data[,c("v", "j")]
                 mat <- reshape2::acast(v~j, data=mat)
                 mat <- mat[order_v, order_j]
                 library(circlize)
                 grid.colors <- 
                   data.frame(names = c(rownames(mat), colnames(mat)),
                              colors= c(colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")))(nrow(mat)), 
                                        colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greens")))(ncol(mat)))
                   )
                 rownames(grid.colors)=grid.colors$names
                 grid.colors$names=NULL
                 grid.colors <- grid.colors %>% as.matrix()
                 
                 
                 circlize::circos.clear()
                 circlize::chordDiagram(mat, grid.col = grid.colors[,"colors"] ,annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
                 circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
                   circlize::circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
                 }, bg.border = NA) 
                 
                 dev.off()
                 
                 
               })
    
  }
  
  if(is.null(arrange.by)){
    # Filter vdj by top and constant region
    if(select.c == "TRB" | select.c == "TRD"){
      data <- 
        data %>% 
        dplyr::filter(c==select.c) %>% 
        dplyr::select(c,v,d,j, barcodes) %>% 
        dplyr::mutate(reAr=paste0(v,"_",d,"_",j))
      
      if(nrow(data)>top){
        
        select.top <- 
          data %>% 
          dplyr::group_by(reAr) %>% 
          dplyr::count() %>% 
          dplyr::arrange(desc(n)) %>% 
          head(top) %>% 
          dplyr::pull(reAr)
        
        data <- 
          data %>% 
          dplyr::filter(reAr %in% select.top)}
      
    }else{
      
      data <- 
        data %>% 
        dplyr::filter(c==select.c) %>% 
        dplyr::select(c,v,j,barcodes) %>% 
        dplyr::mutate(reAr=paste0(v,"_",j))
      
      if(nrow(data)>top){
        
        select.top <- 
          data %>% 
          dplyr::group_by(reAr) %>% 
          dplyr::count() %>% 
          dplyr::arrange(desc(n)) %>% 
          head(top) %>% 
          dplyr::pull(reAr)
        
        data <- 
          data %>% 
          dplyr::filter(reAr %in% select.top)
        
        
      }
      
    }
    
    data %>% head()
    
    
    mat <- data[,c("v", "j")]
    names(mat) <- c("from", "to")
    mat <- reshape2::acast(from~to, data=mat)
    
    
    grid.colors <- 
      data.frame(names = c(rownames(mat), colnames(mat)),
                 colors= c(colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(nrow(mat)), 
                           colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(ncol(mat)))
      )
    rownames(grid.colors)=grid.colors$names
    grid.colors$names=NULL
    grid.colors <- grid.colors %>% as.matrix()
    
    
    circlize::circos.clear()
    circlize::chordDiagram(mat, grid.col = grid.colors[,"colors"] ,annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
    circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
      circlize::circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    }, bg.border = NA) 
    
    
  }
  
  
}




#' @title plotVDJArrangments
#' @author Dieter Henrik Heiland
#' @description SPATA-Immuno
#' @return 
#' @examples 
#' 
#' @export
#' 
plotVDJArrangments<- function(object, 
                              select.c="TRB", 
                              arrange.by=NULL, 
                              top=20, 
                              output=getwd()){
  
  data <- SPATAImmune::GetVDJ(object)
  if(is.null(data)){data <- SPATAImmune::GetSPTCR(object)}
  if(is.null(data)) stop ("No VDJ or SPTCR slot in object")
  
  
  
  
  # If required select spots
  if(!is.null(arrange.by)){
    bc.join <- joinWithFeatures(object, features = arrange.by) %>% select(barcodes, {{arrange.by}})
    data <- 
      data %>% 
      dplyr::left_join(.,bc.join, by="barcodes") %>% 
      as.data.frame()
    
    plot.out <- purrr::map(.x=unique(data[,arrange.by]),
               .f=function(select){
                 data <- data %>% dplyr::filter(!!sym(arrange.by)=={{select}})
                 
                 if(select.c == "TRB" | select.c == "TRD"){
                   data <- 
                     data %>% 
                     dplyr::filter(c==select.c) %>% 
                     dplyr::select(v,d,j, barcodes) %>% 
                     dplyr::mutate(reAr=paste0(v,"_",d,"_",j))
                   
                   if(nrow(data)>top){
                     
                     select.top <- 
                       data %>% 
                       dplyr::group_by(reAr) %>% 
                       dplyr::count() %>% 
                       dplyr::arrange(desc(n)) %>% 
                       head(top) %>% 
                       dplyr::pull(reAr)
                     
                     data <- 
                       data %>% 
                       dplyr::filter(reAr %in% select.top)}
                   
                 }else{
                   
                   data <- 
                     data %>% 
                     dplyr::filter(c==select.c) %>% 
                     dplyr::select(v,j,barcodes) %>% 
                     dplyr::mutate(reAr=paste0(v,"_",j))
                   
                   if(nrow(data)>top){
                     
                     select.top <- 
                       data %>% 
                       dplyr::group_by(reAr) %>% 
                       dplyr::count() %>% 
                       dplyr::arrange(desc(n)) %>% 
                       head(top) %>% 
                       dplyr::pull(reAr)
                     
                     data <- 
                       data %>% 
                       dplyr::filter(reAr %in% select.top)
                     
                     
                   }
                   
                 }
                 
                 plot.out <- SPATAImmune::plotDonut(table(data$reAr))
                 return(plot.out)
               })

    }
  
  if(is.null(arrange.by)){
    # Filter vdj by top and constant region
    if(select.c == "TRB" | select.c == "TRD"){
      data <- 
        data %>% 
        dplyr::filter(c==select.c) %>% 
        dplyr::select(v,d,j, barcodes) %>% 
        dplyr::mutate(reAr=paste0(v,"_",d,"_",j))
      
      if(nrow(data)>top){
        
        select.top <- 
          data %>% 
          dplyr::group_by(reAr) %>% 
          dplyr::count() %>% 
          dplyr::arrange(desc(n)) %>% 
          head(top) %>% 
          dplyr::pull(reAr)
        
        data <- 
          data %>% 
          dplyr::filter(reAr %in% select.top)}
      
    }else{
      
      data <- 
        data %>% 
        dplyr::filter(c==select.c) %>% 
        dplyr::select(v,j,barcodes) %>% 
        dplyr::mutate(reAr=paste0(v,"_",j))
      
      if(nrow(data)>top){
        
        select.top <- 
          data %>% 
          dplyr::group_by(reAr) %>% 
          dplyr::count() %>% 
          dplyr::arrange(desc(n)) %>% 
          head(top) %>% 
          dplyr::pull(reAr)
        
        data <- 
          data %>% 
          dplyr::filter(reAr %in% select.top)
        
        
      }
      
    }
    plot.out <- SPATAImmune::plotDonut(table(data$reAr))
    return(plot.out)
    
    
    
  }
  
  return(plot.out)
  
}


#' @title plotAAlogo
#' @author Dieter Henrik Heiland
#' @description SPATA-Immuno
#' @return 
#' @examples 
#' 
#' @export
#' 
plotAAlogo <- function(object, motif.length=8, cluster.select=2, constant="TRB") {
  
  data <- SPATAImmune::GetVDJ(object, constant=constant)
  
  if(!any(names(data)=="cluster"))stop("No clusters called")
  
  
  logo <- 
    data %>% 
    dplyr::select(CDR3_aa, cluster, AA_length) %>% 
    dplyr::filter(cluster==cluster.select) %>% 
    dplyr::filter(AA_length>motif.length)
  
  #Select only the first 8
  motif.consesus <- map_df(.x=1:length(logo$CDR3_aa), .f=function(i){
    k <- logo$CDR3_aa[i]
    motif.select <- stringr::str_split(k, pattern="") %>% unlist()
    motif.select <- motif.select[1:motif.length]
    out <- data.frame(motif.select) %>% t()
    rownames(out) <- paste0("Seq_",i)
    names(out) <- paste0("AA_",1:motif.length)
    out <- out %>% as.data.frame()
    return(out)
  }) %>% as.data.frame()
  
  flatten.df <- function(df){
    map(.x=names(df), .f=function(tab){df[,tab]}) %>% unlist()
  }
  
  AA <- unique(motif.consesus %>%flatten.df()) 
  mat.motif <- matrix(0, length(AA), motif.length)
  rownames(mat.motif) <- AA
  
  for(i in 1:ncol(mat.motif)){
    
    for(ii in 1:nrow(mat.motif)){
      val <- motif.consesus[i]==rownames(mat.motif)[ii]
      quant <- which(val==T) %>% length()/length(val)
      mat.motif[ii,i] <- quant
    }
  }
  
  #for(i in 1:ncol(mat.motif)){ mat.motif[,i] <- scales::rescale(mat.motif[,i], c(0,1)) }
  
  mat.motif <- reshape2::melt(mat.motif)
  
  require(gglogo)
  p=ggplot(data=mat.motif)+
    gglogo::geom_logo(aes(x=Var2, y=value, label=Var1, fill=Var1),alpha = 0.6)  +
    theme_void()+Seurat::NoLegend()
  
  return(p)
}



#' @title plotBarCompare
#' @author Dieter Henrik Heiland
#' @description SPATA-Immuno
#' @return 
#' @examples 
#' 
#' @export
#' 
plotBarCompare <- function (plot_df, Feat, compare_to, order=NULL) {
  plot_df = plot_df[, c(Feat, compare_to)]
  names(plot_df) = c("Cx", "Tx")
  
  if(!is.null(order)){
    plot_df$Tx <- as.factor(plot_df$Tx)
    p = ggplot(data = plot_df, aes(x = factor(Tx, levels = order), y = 1, fill = Cx)) + 
      geom_bar(position = "fill", stat = "identity") + theme_classic() + 
      xlab(compare_to) + ylab("Percentage") + theme(plot.margin = margin(t = 50, 
                                                                         r = 100, b = 50, l = 100, unit = "pt"), axis.text.y = element_text(color = "black"), 
                                                    axis.text.x = element_text(color = "black", angle = 75, 
                                                                               vjust = 0.5))
  }else{
    
    p = ggplot(data = plot_df, aes(x = Tx, y = 1, fill = Cx)) + 
      geom_bar(position = "fill", stat = "identity") + theme_classic() + 
      xlab(compare_to) + ylab("Percentage") + theme(plot.margin = margin(t = 50, 
                                                                         r = 100, b = 50, l = 100, unit = "pt"), axis.text.y = element_text(color = "black"), 
                                                    axis.text.x = element_text(color = "black", angle = 75, 
                                                                               vjust = 0.5))
  }
  
 
  return(p)
}



#' @title plotAggregatedDimred
#' @author Dieter Henrik Heiland
#' @description SPATA-Immuno
#' @return 
#' @examples 
#' 
#' @export
#'
plotAggregatedDimred <- function(object, factor=20, pt.alpha=.8, pt.size=0.5, return.coords=F,constant){
  #check if UMAP is avaiable
  sample <- object@samples
  data <- SPATAImmune::GetVDJ(object,constant)
  
  if(!any(names(data) %in% c("UMAP_1")))stop(" UMAP not build, please run the runSPTCRdeep function before")
  
  numberofclusters <- length(unique(data$cluster))
  
  # Create random aggregation vectors
  random.x <- runif(numberofclusters, factor,factor*factor) %>% round()
  random.y <- runif(numberofclusters, factor,factor*factor) %>% round()
  
  data.plot <- 
    data_frame(barcodes=data$barcodes, 
               Cluster=data$cluster %>% as.numeric(),
               UMAP_1= data$UMAP_1,
               UMAP_2= data$UMAP_2) %>% 
    mutate(Dim.x = map_chr(.x=.$Cluster, function(c){random.x[as.numeric(c)+1]}) %>% as.numeric(),
           Dim.y = map_chr(.x=.$Cluster, function(c){random.y[as.numeric(c)+1]}) %>% as.numeric()) %>% 
    mutate(Dim1=UMAP_1+jitter(Cluster+Dim.x,3),
           Dim2=UMAP_2+jitter(Cluster+Dim.y,3)) %>% 
    mutate(Cluster=as.character(data$cluster))
  
  
  #Plot data 
  p <-   
    ggplot()+
    geom_point(data=data.plot, mapping=aes(x=Dim1,y=Dim2, color=Cluster), size=pt.size, alpha=pt.alpha)+
    theme_classic()
  
  object@data[[sample]]$VDJ$UMAP_Aggregated <-data.plot 
  
  if(return.coords==T){return(object)}else{return(p)}
  
  
  
}


#' @title plot3DClonality2ST
#' @author Dieter Henrik Heiland
#' @description SPATA-Immuno
#' @return 
#' @examples 
#' 
#' @export
#' 

plot3DClonality2ST <- function(object, color_by, pt.bot=3, pt.up=1, random.spots=10,line.size=0.3, select.cluster=NULL, constant="TRB", angle=list(x=1.87, y=0.88, z=-0.64)){
  
  #get Surface data
  of_sample <- object@samples
  surface <- SPATA2::plotSurface(object, color_by=color_by)$data
  names(surface)[names(surface)==color_by] <- "feat"
  surface$z=0
  
  #Get Aggregated UMAP
  clones <- object@data[[of_sample]]$VDJ$UMAP_Aggregated
  select.constant <- which(SPATAImmune::GetVDJ(object)$c==constant)
  clones <- clones[select.constant, ]
  clones$z=10
  
  
  
  
  # Bring the Cords into the sampe space
  surface$x <- scales::rescale(surface$x, c(0,1))
  surface$y <- scales::rescale(surface$y, c(0,1))
  
  clones$Dim1 <- scales::rescale(clones$Dim1 , c(0,1))
  clones$Dim2 <- scales::rescale(clones$Dim2 , c(0,1))
  
  
  #Layout parameters 
  
  az <- list(
    title ="",
    range = c(0,15),
    showgrid=F,
    zeroline=F,
    showline=F,
    autotick=F,
    ticks= '',
    showticklabels=F
  )
  ax <- list(
    title ="",
    showgrid=F,
    zeroline=F,
    showline=F,
    autotick=F,
    ticks= '',
    showticklabels=F
  )
  ay <- list(
    title ="",
    showgrid=F,
    zeroline=F,
    showline=F,
    autotick=F,
    ticks= '',
    showticklabels=F
  )
  
  
  
  # Create plotly plots
  plot_plotly <-
    plotly::plot_ly(data=clones, 
                    x= ~Dim1, 
                    y= ~Dim2, 
                    z= ~z, 
                    color = ~Cluster,
                    marker=list(#colorscale = "Plasma",
                      reversescale =T,
                      size=pt.up,
                      alpha=0.8),
                    type = "scatter3d") %>% 
    plotly::add_trace(data=surface, 
                      x= ~x, 
                      y= ~y, 
                      z= ~z,
                      color = ~feat,
                      marker=list(#colorscale = "Plasma",
                        reversescale =T,
                        size=pt.bot,
                        alpha=1),
                      type = "scatter3d") %>% 
    plotly::layout(scene = list(zaxis=az, yaxis=ay, xaxis=ax,camera=list(eye = angle)),
                   showlegend = FALSE)
  
  
  #Add positions
  spots <- 
    clones %>% 
    dplyr::left_join(.,surface, by="barcodes")
  
  #Select only randomly 100 points per cluster
  random.select <- map(.x=unique(spots$Cluster), .f=function(c){ 
    all <- which(spots$Cluster == c )
    if(length(all)<random.spots){size=length(all)}else{size=random.spots}
    all %>% sample(., size=size) }) %>% unlist()
  spots <- spots[random.select, ]
  
  
  
  if(!is.null(select.cluster)){
    spots <- 
      clones %>% 
      dplyr::left_join(.,surface, by="barcodes") %>% 
      filter(Cluster %in% {{select.cluster}})
    random.select <- map(.x=unique(spots$Cluster), .f=function(c){ 
      all <- which(spots$Cluster == c )
      if(length(all)<random.spots){size=length(all)}else{size=random.spots}
      all %>% sample(., size=size) }) %>% unlist()
    spots <- spots[random.select, ]
  }
  
  
  
  
  
  x1=spots$Dim1
  x2=spots$x
  y1=spots$Dim2
  y2=spots$y
  z1=spots$z.x
  z2=spots$z.y
  
  #Run the loop to add lines
  
  for(i in 1:length(x1)){
    
    plot_plotly <- plot_plotly %>% 
      plotly::add_paths(x = c(x1[i],x2[i]), 
                        y = c(y1[i],y2[i]), 
                        z = c(z1[i],z2[i]), 
                        color='black', 
                        line = list(width = line.size, alpha=0.1, color = "black", reverscale = FALSE)
      )
    
  }
  
  return(plot_plotly)
  
}

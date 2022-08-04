#' @title runKmers -Deprecated - 
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @return SPATA object
#' @examples 
#' @export
runKmers <- function(seq, kmer.length=4, kmer.nr=NULL, overlap=1){
  
  length.seq <- stringi::stri_length(seq)
  max.seq.length <- max(length.seq, na.rm = T)
  min.seq.length <- min(length.seq, na.rm = T)
  
  if(is.null(kmer.nr)){
    kmer.nr=min.seq.length/4
  }else{kmer.nr=kmer.nr}
  message(paste0(Sys.time(), " --- Sequence range between: ",min.seq.length, " to ", max.seq.length, " aa ---" ))
  
  kmers <- map(.x=seq, .f=function(s){
    #print(s)
    seq.x <- unlist(strsplit(s, "") )
    nr.kmers <- 1+(length(seq.x)/(kmer.length-overlap)) %>% round()
    v=seq(1,(length(seq.x)-kmer.length), by=c(kmer.length-overlap))
    x <- map(.x=v, .f=function(i){seq.x[i:(i+kmer.length)] %>% paste0(collapse="")}) %>% unlist() %>% paste(collapse = " ")
    
    
  }) %>% unlist()
  
}

#' @title modelAA2vec -Deprecated - 
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @return SPATA object
#' @examples 
#' @export
modelAA2vec <- function(){
  
  require(keras)
  
  input_target <- keras::layer_input(shape = 1)
  input_context <- keras::layer_input(shape = 1)
  
  embedding <- keras::layer_embedding(
    input_dim = tokenizer$num_words+1, 
    output_dim = embedding_size, 
    input_length = 1, 
    name = "embedding")
  
  target_vector <- 
    input_target %>% 
    embedding() %>% 
    keras::layer_flatten(name = "Test")
  
  
  context_vector <- 
    input_context %>%
    embedding() %>%
    keras::layer_flatten()
  
  
  dot_product <- keras::layer_dot(list(target_vector, context_vector), axes = 1)
  output <- keras::layer_dense(dot_product, units = 1, activation = "sigmoid")
  
  model <- keras::keras_model(list(input_target, input_context), output)
  model %>% keras::compile(loss = "binary_crossentropy", optimizer = "adam")
  
  return(model)
  
} 



#' @title runSPTCRdeep -Deprecated - 
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @return SPATA object
#' @examples 
#' @export
runSPTCRdeep <- function(object, 
                         constant=NULL,
                         run.AA2vector=T, 
                         run.ANN=T,
                         kmer.length=3,
                         embedding_size = 10,
                         skip_window = 5,
                         epochs.w2v=80,
                         LSTM=T,
                         epochs.LSTM=40,
                         VAE=F,
                         num_sampled = 1,
                         epochs=40, 
                         LSTMBatch=1000,
                         temp.save=T,
                         temp.folder=getwd()){
  
  
  
  
  
  SPATAImmune::verbose.text("The function is Deprecated .... Prepare data")
  if(is.null(constant)){data <- SPATAImmune::GetVDJ(object)}else{data <- SPATAImmune::GetVDJ(object) %>% dplyr::filter(c==constant)}

  
  library(keras)
  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
  K <- keras::backend()
  

  
  #Check if data are avaiable

# Step 1 AA2Vector --------------------------------------------------------
  
  if(run.AA2vector==T){
    
    data <- 
      data %>% 
      dplyr::mutate(AA_length= stringi::stri_length(CDR3_aa)) %>%  
      dplyr::filter(AA_length>kmer.length*2)
    data$CDR3_aa <- stringr::str_replace_all(data$CDR3_aa, pattern="[*]", "X")
    
    kmers <- SPATAImmune::runKmers(data$CDR3_aa, kmer.length=kmer.length, overlap = 1)
    kmers <- base::iconv(kmers)
    
    
    
    
    SPATAImmune::verbose.text("Train AA2vector AA Embedding")
    SPATAImmune::verbose.text("Run Kmers Embedding")
    
    # AA2vector ---------------------------------------------------------------
    
    SPATAImmune::verbose.text("Run Tokenizer ...")
    tokenizer <- 
      keras::text_tokenizer(200000) %>% 
      keras::fit_text_tokenizer(kmers)
    
    tokenizer <- 
      keras::text_tokenizer(tokenizer$word_index %>% length()) %>% fit_text_tokenizer(kmers)
   
    skipgrams_generator <- function(kmers, tokenizer, window_size, negative_samples) {
      
      
      require(reticulate)
      require(purrr)
      gen <- texts_to_sequences_generator(tokenizer, sample(kmers))
      function() {
        skip <- generator_next(gen) %>%
          skipgrams(
            vocabulary_size = tokenizer$num_words, 
            window_size = window_size, 
            negative_samples = 1
          )
        x <- transpose(skip$couples) %>% purrr::map(. %>% unlist %>% as.matrix(ncol = 1))
        y <- skip$labels %>% as.matrix(ncol = 1)
        list(x, y)
      }
    }
    
    
    SPATAImmune::verbose.text("Load Model  ...")

      require(keras)
      
      input_target <- keras::layer_input(shape = 1)
      input_context <- keras::layer_input(shape = 1)
      
      embedding <- keras::layer_embedding(
        input_dim = tokenizer$num_words+1, 
        output_dim = embedding_size, 
        input_length = 1, 
        name = "embedding")
      
      target_vector <- 
        input_target %>% 
        embedding() %>% 
        keras::layer_flatten(name = "Test")
      
      
      context_vector <- 
        input_context %>%
        embedding() %>%
        keras::layer_flatten()
      
      
      dot_product <- keras::layer_dot(list(target_vector, context_vector), axes = 1)
      output <- keras::layer_dense(dot_product, units = 1, activation = "sigmoid")
      
      model <- keras::keras_model(list(input_target, input_context), output)
      model %>% keras::compile(loss = "binary_crossentropy", optimizer = "adam")

    SPATAImmune::verbose.text("Fit Generator  ...")
    model %>% 
      keras::fit_generator(
        skipgrams_generator(kmers, tokenizer, skip_window, negative_samples=1), 
        steps_per_epoch = 10, 
        epochs = epochs.w2v,
        verbose = T)
    
    SPATAImmune::verbose.text("Get Data from AA2verctor...")
    embedding_matrix <- get_weights(model)[[1]]
    words <- 
      data_frame(word = names(tokenizer$word_index), id = as.integer(unlist(tokenizer$word_index))) %>% 
      dplyr::filter(id <= tokenizer$num_words) %>% 
      dplyr::arrange(id)
    
    row.names(embedding_matrix) <- c("UNK", words$word)
    SPATAImmune::verbose.text("Save Data ...")
  #Save data
  if(temp.save==T){
    path=paste0(temp.folder, "/embedding_matrix.RDS")
    saveRDS(embedding_matrix, path)
  }
    
    
  }else{
    SPATAImmune::verbose.text("Load  Embedding matrix ...")
    path=paste0(temp.folder, "/embedding_matrix.RDS")
    embedding_matrix <- readRDS(path)
    }
  
  

    
    
# Step 3 ANN  ---------------------------------------------------------------
    
    
    
    if(run.ANN==T){
    
      if(run.AA2vector==F){
        
        data <- 
          data %>% 
          dplyr::mutate(AA_length= stringi::stri_length(CDR3_aa)) %>%  
          dplyr::filter(AA_length>kmer.length*2)
        data$CDR3_aa <- stringr::str_replace_all(data$CDR3_aa, pattern="[*]", "X")
        
        kmers <- SPATAImmune::runKmers(data$CDR3_aa, kmer.length=kmer.length, overlap = 1)
        kmers <- base::iconv(kmers)
        }
      
      
      
      SPATAImmune::verbose.text("Prepare Data ... ")
      max_len <- map(.x=kmers, .f=function(k){length=str_split(k, pattern=" ") %>% unlist %>% length()}) %>% unlist() %>% max()  

      embedding_seq <- map(.x=kmers, .f=function(k){
        mat <- embedding_matrix[stringr::str_split(k, pattern=" ") %>% unlist() %>%  base::tolower() , ] 
        mat <- rbind(mat, matrix(0, max_len-nrow(mat), ncol(mat)))
      })
      
      array <- SPATAImmune::hlpList2Array(embedding_seq) %>% aperm(., c(3,1,2))
      
      
      if (tensorflow::tf$executing_eagerly())
        tensorflow::tf$compat$v1$disable_eager_execution()
      K <- keras::backend()
      
      SPATAImmune::verbose.text("Run ANN ... ")
      
      
      
      
    if(LSTM==T){
      # LSTM AE model ---------------------------------------------------------------
      
      timesteps <- max_len
      model.LSTM.encoder <- 
        keras::keras_model_sequential() %>% 
        keras::layer_masking(input_shape=c(timesteps, embedding_size), mask_value=0) %>% 
        keras::layer_lstm(units=254, input_shape=c(timesteps, embedding_size), activation="relu", return_sequences = T) %>% 
        keras::layer_lstm(units=128, activation="relu") %>% 
        keras::layer_dense(units = 128, name ="latent_space") %>% 
        keras::layer_repeat_vector(timesteps) %>% 
        keras::layer_lstm(units=timesteps ,return_sequences=T, activation="relu") %>% 
        keras::time_distributed(layer = layer_dense(units = embedding_size))
      
      
      model.LSTM.encoder  %>% keras::compile(loss = 'mse', optimizer = 'adam')
      SPATAImmune::verbose.text("Training LSTM AE ... ")
      model.LSTM.encoder %>% keras::fit(array, array, epochs=epochs.LSTM, batch_size=LSTMBatch)
      
      
      #get latent space
      get.latent.space.2D.LSTM.Autoencoder <- function(model,array){
        keras::keras_model(inputs = model$input, 
                           outputs = keras::get_layer(model, "latent_space")$output) %>% 
          predict(., array)
        
      } 
      
      ls <- get.latent.space.2D.LSTM.Autoencoder(model=model.LSTM.encoder, array)
      #ls %>% dim()
      #umap <- umap::umap(ls)$layout
      #plot(umap, pch=".")
      
      if(temp.save==T){
        path=paste0(temp.folder, "/ls.RDS")
        saveRDS(ls, path)
      }
      
      
    }
    if(VAE==T){
      # VAE model ---------------------------------------------------------------
      
      if (tensorflow::tf$executing_eagerly())
        tensorflow::tf$compat$v1$disable_eager_execution()
      
      library(keras)
      K <- keras::backend()
      
      input.dim=dim(array)
      epsilon_std <- 1.0
      latent_dim <- 32
      
      #Encoder
      input <- keras::layer_input(shape=input.dim[2:3])
      
      encoder <- 
        input %>% 
        keras::layer_flatten() %>% 
        keras::layer_dense(units = 1024, activation = "relu") %>%
        keras::layer_dense(units = 512, activation = "relu") %>%
        keras::layer_dense(units = 256, activation = "relu") %>%
        keras::layer_dense(units = 64, activation = "relu")
      
      z_mean <- layer_dense(encoder, latent_dim)
      z_log_var <- layer_dense(encoder, latent_dim)
      
      sampling <- function(arg){
        z_mean <- arg[, 1:(latent_dim)]
        z_log_var <- arg[, (latent_dim + 1):(2 * latent_dim)]
        
        epsilon <- k_random_normal(
          shape = c(k_shape(z_mean)[[1]]), 
          mean=0.,
          stddev=epsilon_std
        )
        
        z_mean + k_exp(z_log_var/2)*epsilon
      }
      
      z <- layer_concatenate(list(z_mean, z_log_var), name = "latent_space") %>% 
        layer_lambda(sampling)
      
      #Decoder
      
      decoder <- 
        z %>% 
        keras::layer_dense(units = 64, activation = "relu") %>%
        keras::layer_dense(units = 256, activation = "relu") %>%
        keras::layer_dense(units = 512, activation = "relu") %>%
        keras::layer_dense(units = 1024, activation = "relu") %>%
        keras::layer_dense(units = c(input.dim[2]*input.dim[3]), activation = "sigmoid") %>%
        keras::layer_reshape(target_shape=input.dim[2:3])
      
      
      input_size <- c(input.dim[2]*input.dim[3])
      
      vae_loss <- function(input, decoder){
        xent_loss=input_size*loss_binary_crossentropy(input, decoder)
        kl_loss=-0.5*k_mean(1+z_log_var-k_square(z_mean)-k_exp(z_log_var), axis=-1)
        xent_loss + kl_loss
      }
      
      vae <- keras::keras_model(input, decoder)
      vae %>% keras::compile(optimizer = "rmsprop", loss = "binary_crossentropy")
      
      
      vae %>% keras::fit(
        array, array, 
        shuffle = TRUE, 
        epochs = epochs, 
        batch_size = 100, 
      )
      
      
      
      
      #get latent space
      get.latent.space.VAE <- function(model,array){
        keras::keras_model(inputs = model$input, 
                           outputs = keras::get_layer(model, "latent_space")$output) %>% 
          predict(., array)
        
      } 
      
      ls <- get.latent.space.VAE(model=vae, array)
      
      if(temp.save==T){
        path=paste0(temp.folder, "/ls.RDS")
        saveRDS(ls, path)
      }

    }
    
    }else{
      SPATAImmune::verbose.text("Load  Latent Space ...")
      path=paste0(temp.folder, "/ls.RDS")
      ls <- readRDS(path)
    }
    
    
    
    
    
# Step 4 Cluster Analysis --------------------------------------------------------
    SPATAImmune::verbose.text("Run Cluster Analysis")
    #Cluster Analysis
    inp <- ls %>%  t()
    colnames(inp) <- paste0("V_",1:ncol(inp))
    rownames(inp) <- paste0("F_",1:nrow(inp))
    a <- Seurat::CreateSeuratObject(counts = inp)
    a <- Seurat::FindVariableFeatures(a, verbose=F) %>% Seurat::ScaleData(verbose=F)
    a <- Seurat::RunPCA(a, npcs = 50, ndims.print = 1)
    SPATAImmune::verbose.text("Run UMAP")
    a <- Seurat::RunUMAP(a, dims=1:10,verbose=F)
    SPATAImmune::verbose.text("Run Cluster Analysis")
    a <- Seurat::FindNeighbors(a,verbose=F) %>% Seurat::FindClusters(verbose=F)
    
    
    SPATAImmune::verbose.text("Get Data Cluster Analysis")
    # T cell clusters
    data$cluster <- a@meta.data$seurat_clusters
    data$UMAP_1 <- a@reductions$umap@cell.embeddings[,1]
    data$UMAP_2 <- a@reductions$umap@cell.embeddings[,2]
    
    
    #Quantify mean 
    umap.mean <- 
      data %>% 
      dplyr::select(cluster,UMAP_1, UMAP_2) %>% 
      dplyr::group_by(cluster) %>% 
      dplyr::summarise_all(mean) %>% 
      mutate(size=data %>% dplyr::count(cluster) %>% dplyr::pull(n))
    
    #plot.umap.sum <- ggplot(data=umap.mean, aes(x=UMAP_1, y=UMAP_2, size=size, color=cluster))+geom_point()+theme_classic()+Seurat::NoLegend()
    
    
    if(is.null(constant)){
      
      #Create motives of each cluster
      SPATAImmune::verbose.text("Export Data Cluster Analysis")
      sample <- object@samples
      object@data[[sample]]$VDJ$IgBlast <- data
      object@data[[sample]]$VDJ$UMAP <- umap.mean
      
    }else{
      
      #Create motives of each cluster
      SPATAImmune::verbose.text("Export Data Cluster Analysis")
      sample <- object@samples
      object@data[[sample]]$VDJ$Clonality$TRA <- NULL
      object@data[[sample]]$VDJ$Clonality$TRB <- NULL
      object@data[[sample]]$VDJ$Clonality$TRG <- NULL
      object@data[[sample]]$VDJ$Clonality$TRD <- NULL
    
      object@data[[sample]]$VDJ$Clonality[[constant]]$IgBlast <- data
      object@data[[sample]]$VDJ$Clonality[[constant]]$UMAP <- umap.mean
    
}
    
    return(object)    
  
  
  
  
  
  
  
  
}



#' @title runClonalityAnalysis -Deprecated - 
#' @author Dieter Henrik Heiland
#' @description This function import the VDJ data from the SPTCR-Pipline or IgBlast
#' @return SPATA object
#' @examples 
#' @export
runClonalityAnalysis <- function(object, motif.length=8, constant="TRB"){
  
  data <- SPATAImmune::GetVDJ(object, constant = constant)
  
  if(!any(names(data)=="cluster"))stop("No clusters called")
  
  
  
  # 1. Calculate the Motif Entropy per Cluster ------------------------------
  SPATAImmune::verbose.text(" Run Cluster Validation ")
  val.cluster <- map_df(.x=unique(data$cluster), .f=function(cluster.select){
    
    logo <- 
      data %>% 
      dplyr::select(CDR3_aa, cluster, AA_length) %>% 
      dplyr::filter(cluster==cluster.select)
    
    #Consensus Dataframe
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
    
    AA <- unique(motif.consesus %>% SPATAImmune::hlpFlattenDF()) 
    mat.motif <- matrix(0, length(AA), motif.length)
    rownames(mat.motif) <- AA
    
    for(i in 1:ncol(mat.motif)){
      
      for(ii in 1:nrow(mat.motif)){
        val <- motif.consesus[i]==rownames(mat.motif)[ii]
        quant <- which(val==T) %>% length()/length(val)
        mat.motif[ii,i] <- quant
      }
    }
    
    #remove seq with long repetetive sequences 
    Rep.motif <- c(apply(mat.motif, 1, max) %>% sum())/motif.length 
    
    variance.motif <- c(apply(mat.motif, 2, max) %>% sum())/motif.length
    
    variance.motif <- variance.motif+Rep.motif
    
    #library(gglogo)
    #ggplot(data=reshape2::melt(mat.motif))+
    #  gglogo::geom_logo(aes(x=Var2, y=value, label=Var1, fill=Var1),alpha = 0.6)  +
    #  theme_void()+Seurat::NoLegend()
    
    ret.def <- 
      data.frame(Cluster=as.character(cluster.select),
                 Motif.var.score=variance.motif)
    
    return(ret.def)
    
  })  
  bc_cluster <- map_df(.x=unique(data$cluster), .f=function(cluster.select){
    
    bc_cluster <- 
      data %>% 
      dplyr::select(CDR3_aa, cluster, AA_length, barcodes) %>% 
      dplyr::filter(cluster==cluster.select) %>% 
      dplyr::pull(barcodes) %>% 
      unique() %>% 
      length()
    
    ret.def <- 
      data.frame(Cluster=as.character(cluster.select),
                 bc_cluster=bc_cluster)
    
    return(ret.def)
    
    
  })
  clonality_cluster <- map_df(.x=unique(data$cluster), .f=function(cluster.select){
    
    Arr_cluster <- 
      data %>% 
      dplyr::mutate(Arr = paste0(c,"_", v,"_",d,"_",j)) %>% 
      dplyr::select(CDR3_aa, cluster, AA_length, Arr) %>% 
      dplyr::filter(cluster==cluster.select) %>% 
      dplyr::pull(Arr) %>% 
      unique() %>% 
      length()
    
    ret.def <- 
      data.frame(Cluster=as.character(cluster.select),
                 Arr_cluster=Arr_cluster)
    
    return(ret.def)
    
    
  })
  #Merge
  val.cluster <- 
    val.cluster %>% 
    dplyr::left_join(., bc_cluster, by="Cluster") %>% 
    dplyr::left_join(., clonality_cluster, by="Cluster")
  
  
  
  
  plot.df <- 
    val.cluster %>% 
    arrange(Motif.var.score) %>% filter(Motif.var.score>0.5 ) %>% 
    dplyr:: mutate(spatial.score= bc_cluster %>% scales::rescale(., c(min(Motif.var.score), max(Motif.var.score))) ) %>% 
    dplyr:: mutate(Arr_cluster.score= Arr_cluster %>% scales::rescale(., c(min(Motif.var.score), max(Motif.var.score))) )
  
  
  validation <- list(val.cluster,plot.df)
  names(validation) <- c("Cluster_Validation", "ggplot.tab")
  sample <- object@samples
  object@data[[sample]]$VDJ$ClusterValidation$TRA <- NULL
  object@data[[sample]]$VDJ$ClusterValidation$TRB <- NULL
  object@data[[sample]]$VDJ$ClusterValidation$TRG <- NULL
  object@data[[sample]]$VDJ$ClusterValidation$TRD <- NULL
  
  object@data[[sample]]$VDJ$ClusterValidation[[constant]] <- validation
  
  return(object)
  
  
  
}



#' @title runSpaceXR
#' @author Dieter Henrik Heiland
#' @param object SPATA object.
#' @param cluster Meta data parameter.
#' @param Seurat Seurat data (scRNA-seq).
#' @description Single Cell integration
#' @return SPATA object with integrated clusters (from Seurat)
#' @examples 
#' @export
#'
#'
runSpaceXR <- function(Seurat,cluster, SPATA){
  
  
  library(spacexr)
  library(Matrix)
  
  counts <- Seurat@assays$RNA@counts %>% as.matrix()
  counts[1:10, 1:10]
  meta_data <-  Seurat@meta.data %>% as.data.frame() %>% dplyr::select({{cluster}},nFeature_RNA) # load in meta_data (barcodes, clusters, and nUMI)
  
  cell_types <- meta_data[,cluster]; names(cell_types) <- rownames(meta_data) # create cell_types named list
  
  cell_types <- as.factor(cell_types) # convert to factor data type
  nUMI <- colSums(counts)
  
  ### Create the Reference object
  reference <- Reference(counts, cell_types, nUMI)
  
  
  
  counts <- SPATA2::getCountMatrix(SPATA) %>% as.matrix()
  coords <- SPATA2::getCoordsDf(SPATA) %>% dplyr::select(barcodes, x, y) %>% as.data.frame()
  rownames(coords) <- coords$barcodes
  coords <- coords[,c(2:3)]
  names(coords) <- c("xcoord","ycoord")
  
  nUMI <- colSums(counts) 
  puck <- SpatialRNA(coords, counts, nUMI)
  
  
  myRCTD <- create.RCTD(puck, reference, max_cores = 5)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
  results <- myRCTD@results
  
  norm_weights = normalize_weights(results$weights) 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA
  
  norm_weights <- 
    norm_weights %>% 
    as.data.frame() %>% 
    rownames_to_column("barcodes") %>% 
    left_join(SPATA %>% getCoordsDf() %>% dplyr::select(barcodes),., by="barcodes")
  
  norm_weights[is.na(norm_weights)] <- 0
  
  SPATA <- SPATA %>% SPATA2::addFeatures(.,norm_weights, overwrite = T)
  
  return(SPATA)
  
}

















# Functions for IgBlust import and Reconstruction for VAE -----------------

library(tidyverse)
ReadIgBlast <- function(IgBlast, Ampbinner, pos){
  
  message(paste0(Sys.time()," --- Read in files ---- "))
  
  blast <- read.table(IgBlast, sep="\t", header=T)
  bc <- read.csv(Ampbinner, sep="\t", header=F)

  #filter data and align reads
  bc.corrected <- 
    bc %>% 
    dplyr::select(V1,V2) %>% 
    rename("sequence_id":=V1) %>% 
    rename("bc":=V2)
  
  blast.test=blast %>% head()
  
  blast.filter <- 
    blast %>% 
    dplyr::filter(complete_vdj==T) %>%
    #tidyr::separate(.,col=sequence_id, into=c("X", "sequence_id"), sep="[|]") %>% 
    #tidyr::separate(.,col=sequence_id, into=c("sequence_id", "X"), sep="strand") %>% 
    #dplyr::select(-X) %>% 
    dplyr::left_join(., bc.corrected, by="sequence_id")
  blast.filter.2 <- blast.filter %>% dplyr::filter(!is.na(cdr2))
  
  
  
  
}
list2ary = function(input.list){  #input a list of lists
  rows.cols <- dim(input.list[[1]])
  sheets <- length(input.list)
  output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
  colnames(output.ary) <- colnames(input.list[[1]])
  row.names(output.ary) <- row.names(input.list[[1]])
  return(output.ary)    # output as a 3-D array
}
CreateShape <- function(blast.out, var="cdr2", filter=NULL, align=NULL){
  #The data shape input
  #CD3
  length.r <- stringi::stri_length(blast.out[,var])
  max.nuc.length <- max(stringi::stri_length(blast.out[,var]), na.rm = T)
  
  if(!is.null(filter)){
    blast.out <- blast.out[which(length.r>filter), ]
  }
  
  if(!is.null(align)){
    
    blast.out <- blast.out %>% dplyr::left_join(data.frame(sequence_id=align), . , by="sequence_id")

  }
  
  
  seq_ID <- blast.out$sequence_id
  
  
  #synthezide CDR3 with similar length into an array
  #G		Glycine		Gly									P		Proline		Pro
  #A		Alanine		Ala									V		Valine		Val
  #L		Leucine		Leu									I		Isoleucine		Ile
  #M		Methionine		Met									C		Cysteine		Cys
  #F		Phenylalanine		Phe									Y		Tyrosine		Tyr
  #W		Tryptophan		Trp									H		Histidine		His
  #K		Lysine		Lys									R		Arginine		Arg
  #Q		Glutamine		Gln									N		Asparagine		Asn
  #E		Glutamic Acid		Glu									D		Aspartic Acid		Asp
  #S		Serine		Ser									T		Threonine		Thr
  
  nuc.vector <- c("G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T")
  CDR3mat <- matrix(0, max.nuc.length, 20)
  colnames(CDR3mat) <- nuc.vector
  
  #CDR3 to mat
  #parralel processing
  
  options(future.globals.maxSize=20000*1024^2)
  future::plan("multiprocess", workers=8)
  
  message(paste0(Sys.time()," --- Process Shape data :",var ," ---- "))
  
  list.cdr3 <- 
    furrr::future_map(.x=blast.out[,var], .f=function(r){
      
      if(is.na(r)){
        
        r.info <- 
          data.frame(nuc=rep(0,max.nuc.length),
                     nr=1:max.nuc.length)
        
      }else{
        
        
        nuc <- unlist(strsplit(r, ""))
        nuc <- nuc[!nuc=="*"]
        split <- c(length(nuc)/2) %>% round()
        l <- length(nuc)
        
        check <- c(identical(nuc, character(0)), is.na(nuc))
        
        if(any(check==T)){
          r.info <- 
            data.frame(nuc=rep(0,max.nuc.length),
                       nr=1:max.nuc.length)
        }else{
          if(length(nuc)==1){
            
            r.info <- 
              data.frame(nuc=c(nuc,rep(0, c(max.nuc.length-(l)))),
                         nr=1:max.nuc.length)
            
          }else{
            r.info <- 
              data.frame(nuc=c(
                nuc[1:split], 
                rep(0, c(max.nuc.length-(l))),
                nuc[(split+1):length(nuc)]),
                nr=1:max.nuc.length)
          }
        }
        
        
        
        
      }
      
      
      
      data <-  reshape2::melt(CDR3mat) %>% 
        mutate(value=map2(.x=Var1, .y=c(Var2 %>% as.character()), .f=function(x,y){
          if(any(r.info$nr==x) ==T){
            r.info.sub <- r.info %>% dplyr::filter(nr==x)
            if(any(r.info.sub$nuc==y)==T){out=1}else{(out=0)}
          }else{(out=0)}
          return(out)
        }) %>% unlist()) %>% 
        reshape2::dcast(Var1~Var2, value.var = "value") %>% 
        dplyr::select(-Var1)
      
      
      
    }, .progress=T)
  
  message(paste0(Sys.time()," --- Transforme into array ---"))
  
  array <- list2ary(list.cdr3) %>% aperm(., c(3,1,2))
  
  output <- list(array, seq_ID)
  names(output)=c("array", "sequence_id")
  return(output)
  
  
  
}
#Transform data to input shape
TransformeVAEInput <- function(blast.out){
  
  
  CD3 <- CreateShape(blast.out, var="cdr2", filter=3)
  
  V <- CreateShape(blast.out, var="v_germline_alignment", filter=NULL, align = CD3$sequence_id )
  D <- CreateShape(blast.out, var="d_sequence_alignment", filter=NULL, align = CD3$sequence_id )
  J <- CreateShape(blast.out, var="j_germline_alignment", filter=NULL, align = CD3$sequence_id )
  
  
  ##Merge all channels
  out <- list(CD3$array,V$array, D$array,J$array, CD3$sequence_id )
  names(out)=c("CD3","V", "D","J", "sequence_id")
  return(out)
  
  
  
}

#TEst function

Ampbinner <- "AmpBinner_Big_Test_20_09.all_reads.txt"
IgBlast <- "middle 11.49.59.tsv" 
library(tidyverse)
input.list <- ReadIgBlast(IgBlast,Ampbinner)

blast.out <- input.list %>% head(1000)
blast.test <- input.list %>% tail(1000)

data.use  <- TransformeVAEInput(blast.out)



# VAE Architecture --------------------------------------------------------

ModelTcellVAE <- function(){
  
  # Parameter -------------------------------------------------
  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
  
  library(keras)
  K <- keras::backend()
  epsilon_std <- 1.0

  # Model -------------------------------------------------------------------
  activation.1="tanh"
  activation.2="sigmoid"
  activation.3="relu"
  
  loss.imp <- "kld"
  
  loss.all <- list(CD3Rx = loss.imp, 
                   Vx = loss.imp, 
                   Dx = loss.imp, 
                   Jx = loss.imp)
  
  
  #Define Functions
  normalize.layer <- function(tensor){
    out.new <- 
      layer_flatten(tensor) %>% 
      #layer_dense(units = c(shape[1]*shape[2])) %>% 
      layer_dense(units = c(100*100), activation=activation.1) %>% 
      layer_reshape(target_shape = c(100, 100))
    
    return(out.new)
  }
  encoder.run <- function(x, filter, kernel_size){
    out <- 
      x %>% 
      layer_conv_1d(filter, kernel_size=kernel_size, activation=activation.3, padding="same") %>% 
      layer_batch_normalization()
    return(out)
  }
  decoder.run <- function(x, filter, kernel_size){
    out <- 
      x %>% 
      layer_conv_1d(filter, kernel_size=kernel_size, activation=activation.3, padding="same") %>% 
      layer_batch_normalization()
      return(out)
  }
  upscale <- function(tensor, shape, name){
    tensor %>% 
      layer_flatten() %>% 
      layer_dense(units = c(shape[1]*shape[2]),activation=activation.1) %>% 
      layer_reshape(target_shape = shape, name = name)
    
  }
  
  
  
  #create flatten input with all information
  CDR.dim=data.use$CD3 %>% dim()
  V.dim=data.use$V %>% dim()
  D.dim=data.use$D %>% dim()
  J.dim=data.use$J %>% dim()
  
  c.inp <- layer_input(shape = c(CDR.dim[2], CDR.dim[3]))
  v.inp <- layer_input(shape = c(V.dim[2], V.dim[3]))
  d.inp <- layer_input(shape = c(D.dim[2], D.dim[3]))
  j.inp <- layer_input(shape = c(J.dim[2], J.dim[3]))
  
  C <- normalize.layer(c.inp)
  V <- normalize.layer(v.inp)
  D <- normalize.layer(d.inp)
  J <- normalize.layer(j.inp)

  #Encoder
  #input <- layer_input(shape=c(100, 400))
  
  input <- layer_concatenate(list(C, V, D, J))
  
  d1 <- encoder.run(input, 1024, kernel_size=c(5))
  d3 <- encoder.run(d1, 512, kernel_size=c(5)) 
  flat <- layer_flatten(d3)
  #Create latent Space
  dense1 <- layer_dense(flat, units = 256, activation=activation.1)
  
  z_mean <- layer_dense(dense1, units = 16, name="z_mean")
  z_log_var <- layer_dense(dense1, units = 16, name="z_log_var")
  
  sampling <- function(arg){
    z_mean <- arg[, 1:(16)]
    z_log_var <- arg[, (16 + 1):(2 * 16)]
    
    epsilon <- k_random_normal(
      shape = c(k_shape(z_mean)[[1]]), 
      mean=0.,
      stddev=epsilon_std
    )
    
    z_mean + k_exp(z_log_var/2)*epsilon
  }
  
  z <- 
    layer_concatenate(list(z_mean, z_log_var), name = "LS") %>% 
    layer_lambda(sampling)
  
  up.dense3 <- layer_dense(z, units = c(100*256),activation=activation.3)
  reshape.up <- layer_reshape(up.dense3, target_shape = c(100, 256))
  u1 <- decoder.run(reshape.up, 256, kernel_size=c(5))
  u2 <- decoder.run(u1, 512, kernel_size=c(5))
  dec_output <- layer_conv_1d(u2, filter=400, kernel_size=c(3), activation=activation.2, padding="same")

  CD3Rx <- upscale(tensor=dec_output, shape = c(CDR.dim[2], CDR.dim[3]), name="CD3Rx")
  Vx <- upscale(tensor=dec_output,shape = c(V.dim[2], V.dim[3]), name="Vx")
  Dx <- upscale(tensor=dec_output,shape = c(D.dim[2], D.dim[3]), name="Dx")
  Jx <- upscale(tensor=dec_output,shape = c(J.dim[2], J.dim[3]), name="Jx")
  
  

  # Compile -----------------------------------------------------------------

  
  # end-to-end autoencoder
  vae <- keras_model(inputs = c(c.inp,v.inp,d.inp,j.inp), 
                     outputs= c(CD3Rx,Vx,Dx,Jx))
  
  vae %>% compile(optimizer=keras::optimizer_adam(), 
                  loss = loss.all,
                  loss_weights = list(CD3Rx = 1.0, Vx = 0.5, Dx = 0.5, Jx = 0.1)
                  )  
  
  #vae <- keras_model(input, dec_output)
  
  # encoder, from inputs to latent space
  #encoder <- keras_model(c(CD3R,V,D,J), z_mean)
  
  



  
  return(vae)
  
  
  
  
}
ModelTcellVAE.single <- function(){
  
  # Parameter -------------------------------------------------
  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
  
  library(keras)
  K <- keras::backend()
  epsilon_std <- 1.0
  
  # Model -------------------------------------------------------------------
  
  activation.1="tanh"
  
  #Define Functions
  normalize.layer <- function(tensor){
    out.new <- 
      layer_flatten(tensor) %>% 
      #layer_dense(units = c(shape[1]*shape[2])) %>% 
      layer_dense(units = c(100*100), activation=activation.1, kernel_regularizer=regularizer_l1_l2(l1 = 0.01, l2 = 0.01)) %>% 
      layer_reshape(target_shape = c(100, 100))
    
    return(out.new)
  }
  encoder.run <- function(x, filter, kernel_size){
    out <- 
      x %>% 
      layer_conv_1d(filter, kernel_size=kernel_size, activation=activation.1, padding="same")
    return(out)
  }
  decoder.run <- function(x, filter, kernel_size){
    out <- 
      x %>% 
      layer_conv_1d(filter, kernel_size=kernel_size, activation=activation.1, padding="same") %>% 
      return(out)
  }
  

  #create flatten input with all information
  CDR.dim=data.use$CD3 %>% dim()
  c.inp <- layer_input(shape = c(CDR.dim[2], CDR.dim[3]))

  d1 <- encoder.run(c.inp, 512, kernel_size=c(5))
  d2 <- encoder.run(d1, 128, kernel_size=c(5))
  flat <- layer_flatten(d2)
  #Create latent Space
  dense1 <- layer_dense(flat, units = 516, activation = activation.1)
  dense2 <- layer_dense(dense1, units = 256, activation = activation.1)
  z_mean <- layer_dense(dense2, units = 16, name="z_mean")
  z_log_var <- layer_dense(dense2, units = 16, name="z_log_var")
  sampling <- function(arg){
    z_mean <- arg[, 1:(16)]
    z_log_var <- arg[, (16 + 1):(2 * 16)]
    
    epsilon <- k_random_normal(
      shape = c(k_shape(z_mean)[[1]]), 
      mean=0.,
      stddev=epsilon_std
    )
    
    z_mean + k_exp(z_log_var/2)*epsilon
  }
  z <- 
    layer_concatenate(list(z_mean, z_log_var), name = "LS") %>% 
    layer_lambda(sampling)

  up.dense1 <- layer_dense(z, units = c(100*256),activation = activation.1)
  reshape <- layer_reshape(up.dense1, target_shape = c(100, 256))
  u1 <- decoder.run(reshape, 256, kernel_size=c(5))
  dec_output <- layer_conv_1d(u1, filter=400, kernel_size=c(3), activation=activation.1, padding="same")
  
  #return in format of input
  
  upscale <- function(tensor, shape, name){
    tensor %>% 
      layer_flatten() %>% 
      layer_dense(units = c(shape[1]*shape[2]),activation = activation.1) %>% 
      layer_reshape(target_shape = shape, name = name)
    
  }
  CD3Rx <- upscale(tensor=dec_output, shape = c(CDR.dim[2], CDR.dim[3]), name="CD3Rx")

  
  
  
  # Compile -----------------------------------------------------------------
  
  
  # end-to-end autoencoder
  vae <- keras_model(inputs = c(c.inp), 
                     outputs= c(CD3Rx))
  
  vae %>% compile(optimizer=keras::optimizer_adam(clipnorm =1), 
                  loss = 'mse' )  
  
  #vae <- keras_model(input, dec_output)
  
  # encoder, from inputs to latent space
  #encoder <- keras_model(c(CD3R,V,D,J), z_mean)
  
  
  
  
  
  
  return(vae)
  
  
  
  
}
ModelTcellVAE.LS <- function(){
  
  # Parameter -------------------------------------------------
  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
  
  library(keras)
  K <- keras::backend()
  epsilon_std <- 1.0
  
  # Model -------------------------------------------------------------------
  
  
  #Define Functions
  normalize.layer <- function(tensor){
    out.new <- 
      layer_flatten(tensor) %>% 
      #layer_dense(units = c(shape[1]*shape[2])) %>% 
      layer_dense(units = c(100*100), activation="relu", kernel_regularizer=regularizer_l1_l2(l1 = 0.01, l2 = 0.01)) %>% 
      layer_reshape(target_shape = c(100, 100))
    
    return(out.new)
  }
  encoder.run <- function(x, filter, kernel_size){
    out <- 
      x %>% 
      layer_conv_1d(filter, kernel_size=kernel_size, activation="relu", padding="same")
    return(out)
  }
  decoder.run <- function(x, filter, kernel_size){
    out <- 
      x %>% 
      layer_conv_1d(filter, kernel_size=kernel_size, activation="relu", padding="same") %>% 
      return(out)
  }
  
  
  
  
  #create flatten input with all information
  CDR.dim=data.use$CD3 %>% dim()
  V.dim=data.use$V %>% dim()
  D.dim=data.use$D %>% dim()
  J.dim=data.use$J %>% dim()
  
  c.inp <- layer_input(shape = c(CDR.dim[2], CDR.dim[3]))
  v.inp <- layer_input(shape = c(V.dim[2], V.dim[3]))
  d.inp <- layer_input(shape = c(D.dim[2], D.dim[3]))
  j.inp <- layer_input(shape = c(J.dim[2], J.dim[3]))
  
  C <- normalize.layer(c.inp)
  V <- normalize.layer(v.inp)
  D <- normalize.layer(d.inp)
  J <- normalize.layer(j.inp)
  
  #Encoder
  #input <- layer_input(shape=c(100, 400))
  
  input <- layer_concatenate(list(C, V, D, J))
  
  d1 <- encoder.run(input, 512, kernel_size=c(5))
  d2 <- encoder.run(d1, 1024, kernel_size=c(5))
  d3 <- encoder.run(d2, 512, kernel_size=c(5)) 
  d4 <- encoder.run(d3, 256, kernel_size=c(5))
  
  enc_output <- d4
  
  # Add flatten vector into latent space
  
  flat <- layer_flatten(enc_output)
  #Create latent Space
  dense1 <- layer_dense(flat, units = 516, activation = "relu")
  dense2 <- layer_dense(dense1, units = 256, activation = "relu")
  #latent.space <- layer_dense(dense2, units = 16, activation = "relu", name="latentspace")
  
  
  
  
  z_mean <- layer_dense(dense2, units = 16, name="z_mean")
  z_log_var <- layer_dense(dense2, units = 16, name="z_log_var")
  
  sampling <- function(arg){
    z_mean <- arg[, 1:(16)]
    z_log_var <- arg[, (16 + 1):(2 * 16)]
    
    epsilon <- k_random_normal(
      shape = c(k_shape(z_mean)[[1]]), 
      mean=0.,
      stddev=epsilon_std
    )
    
    z_mean + k_exp(z_log_var/2)*epsilon
  }
  
  z <- 
    layer_concatenate(list(z_mean, z_log_var), name = "LS") %>% 
    layer_lambda(sampling)
  
  
  # Reshape Latent space
  
  #up.dense.input <- layer_input(shape = c(2))
  up.dense1 <- layer_dense(z, units = c(256))
  up.dense2 <- layer_dense(up.dense1, units = c(516))
  up.dense3 <- layer_dense(up.dense2, units = c(100*256))
  
  reshape <- layer_reshape(up.dense3, target_shape = c(100, 256))
  
  
  u1 <- decoder.run(reshape, 256, kernel_size=c(5))
  u2 <- decoder.run(u1, 512, kernel_size=c(5))
  u3 <- decoder.run(u2, 1024, kernel_size=c(5))
  u4 <- decoder.run(u3, 512, kernel_size=c(5))
  dec_output <- layer_conv_1d(u4, filter=400, kernel_size=c(3), activation="sigmoid", padding="same")
  
  #return in format of input
  
  upscale <- function(tensor, shape, name){
    tensor %>% 
      layer_flatten() %>% 
      layer_dense(units = c(shape[1]*shape[2])) %>% 
      layer_reshape(target_shape = shape, name = name)
    
  }
  
  
  CD3Rx <- upscale(tensor=dec_output, shape = c(CDR.dim[2], CDR.dim[3]), name="CD3Rx")
  Vx <- upscale(tensor=dec_output,shape = c(V.dim[2], V.dim[3]), name="Vx")
  Dx <- upscale(tensor=dec_output,shape = c(D.dim[2], D.dim[3]), name="Dx")
  Jx <- upscale(tensor=dec_output,shape = c(J.dim[2], J.dim[3]), name="Jx")
  
  
  
  # Compile -----------------------------------------------------------------
  
  
  # end-to-end autoencoder
  vae <- keras_model(inputs = c(c.inp,v.inp,d.inp,j.inp), 
                     outputs= z)
  
  vae %>% compile(optimizer=keras::optimizer_adam(clipnorm =1), 
                  loss = "mse"
  )  
  
  #vae <- keras_model(input, dec_output)
  
  # encoder, from inputs to latent space
  #encoder <- keras_model(c(CD3R,V,D,J), z_mean)
  
  
  
  
  
  
  return(vae)
  
  
  
  
}

# Training ----------------------------------------------------------------
library(keras)
model <- ModelTcellVAE()

train <- list(data.use$CD3[1:400,,], data.use$V[1:400,,], data.use$D[1:400,,], data.use$J[1:400,,])
test <- list(data.use$CD3[401:600,,], data.use$V[401:600,,], data.use$D[401:600,,], data.use$J[401:600,,])


model %>% fit(
  train, 
  train, 
  shuffle = T, 
  epochs = 10,
  batch_size =25,
  validation_data = list(test, test)
)



###### Model for only CDR3 

model <- ModelTcellVAE.single()

train <- list(data.use$CD3)
test <- list(data.use$CD3[401:600,,])


model %>% fit(
  train, 
  train, 
  shuffle = T, 
  epochs = 10,
  batch_size =25,
  validation_data = list(test, test)
)


#model %>% save_model_weights_tf("Test1.ckpt")
#model %>% load_model_weights_tf("Test1.ckpt")


# Reconstruction ----------------------------------------------------------

get.latent.space <- function(model,array){
  keras_model(inputs = model$input, 
              outputs = get_layer(model, "LS")$output) %>% 
    predict(., array)
  
} 

LS <- get.latent.space(model=model, array=list(data.use$CD3, data.use$V, data.use$D, data.use$J))



LS <- get.latent.space(array=list(data.use$CD3))
plot(umap::umap(LS)$layout)


data.pred <- ModelTcellVAE() %>% predict(train)
umap <- get.latent.space(model=model, array=train) %>% umap::umap()
umap$layout %>% as.data.frame() %>% ggplot(data=., aes(x=V1, y=V2))+geom_point(size=1, alpha=0.2)+theme_classic()




#Run TCR Clustering and relative 








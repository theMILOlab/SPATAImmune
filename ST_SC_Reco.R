



# Model for sc st RNA seq transformation ----------------------------------



#Read and create scRNAseq and stRNA-seq datasets


#load scRNA-seq dataset
library(tidyverse)
library(Seurat)
scdata <- readRDS("object.NATNC.RDS")

cell.types <- scdata@meta.data$cluster %>% unique()
most.var.genes <- Seurat::VariableFeatures(scdata)


simulate.STspots=function(scdata,
                          n.spots=1000,
                          cells.per.spot=5,
                          var.genes=most.var.genes)
  {
  

# Prepare Spots -----------------------------------------------------------
  
  simulated.mat <- matrix(0, length(var.genes), n.spots)
  colnames(simulated.mat) <- paste0("ST_", 1:n.spots)
  rownames(simulated.mat) <- var.genes
  
  #Create a random distrubution of cells 
  
  random.cells <- scdata@meta.data %>% rownames() %>% sample(., size=n.spots*cells.per.spot)
  seq.v <- seq(1,c(length(random.cells)-(cells.per.spot-1)), by=cells.per.spot) %>% sample()

  #Create random spots
  for(i in 1:n.spots){
    seq <- seq.v[i]
    #print(i)
    #get cell Count matrix
    mat <- Seurat::GetAssayData(scdata, slot="counts")[var.genes, random.cells[seq:(seq+(cells.per.spot-1))]] %>% as.matrix() %>% t() %>% colMeans()
    simulated.mat[names(mat),i] <- mat
  }
  # Prepare sc Data  ----------------------------------------------------------- 
  #Create gene expression scDataset
  simulate.sc <- map(.x=1:n.spots, .f=function(i){
    seq <- seq.v[i]
    Seurat::GetAssayData(scdata, slot="counts")[var.genes, random.cells[seq:(seq+(cells.per.spot-1))]] %>% as.matrix()
  })
  
  list2ary = function(input.list){  #input a list of lists
    rows.cols <- dim(input.list[[1]])
    sheets <- length(input.list)
    output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
    colnames(output.ary) <- colnames(input.list[[1]])
    row.names(output.ary) <- row.names(input.list[[1]])
    return(output.ary)    # output as a 3-D array
  }
  simulate.sc.array <- list2ary(simulate.sc) %>%  aperm(., c(3,1,2))
  
  out <- list(simulated.mat, simulate.sc.array)
  names(out) <- c("st", "sc")
  
  return(out)
  
  
}


data <- simulate.STspots(scdata, var.genes=most.var.genes)

# Model -------------------------------------------------------------------

if (tensorflow::tf$executing_eagerly())
  tensorflow::tf$compat$v1$disable_eager_execution()

library(keras)
K <- keras::backend()
epsilon_std <- 1.0
n.spots=1000
activation.1="tanh"
activation.2="sigmoid"
activation.3="relu"


#Loss Functions


# Loss function for the VAE
# Loss function comprised of two parts, Cross_entropy, and
# divergence

vae_loss <- function(inputs, finalLayer){
  K <- backend()
  reconstruction_loss = k_sum(k_square(finalLayer - inputs))
  kl_loss = - 0.5 * k_sum(1 + z_log_sigmaFull - k_square(z_meanFull) - k_square(k_exp(z_log_sigmaFull)), axis=-1)
  total_loss = k_mean(reconstruction_loss + kl_loss)
  return(total_loss)
}


# Loss function for the left VAE
# Loss function comprised of two parts, Cross_entropy, and
# divergence
left_vae_loss <- function(inputs, finalLayer){
  K <- backend()
  reconstruction_loss = k_sum(k_square(finalLayer - inputs))
  kl_loss = - 0.5 * k_sum(1 + z_log_sigmaLeft - k_square(z_meanLeft) - k_square(k_exp(z_log_sigmaLeft)), axis=-1)
  total_loss = k_mean(reconstruction_loss + kl_loss)
  return(total_loss)
}


# Loss function for the right VAE
# Loss function comprised of two parts, Cross_entropy, and
# divergence
right_vae_loss <- function(inputs, finalLayer){
  K <- backend()
  reconstruction_loss = k_sum(k_square(finalLayer - inputs))
  kl_loss = - 0.5 * k_sum(1 + z_log_sigmaRight - k_square(z_meanRight) - k_square(k_exp(z_log_sigmaRight)), axis=-1)
  total_loss = k_mean(reconstruction_loss + kl_loss)
  return(total_loss)
}
 

# Left Encoder (scData)
library(keras)
leftEncoderInput <-  layer_input(shape=c(3000, 5), dtype="float32")
leftEncoderLayer1 <- layer_conv_1d(leftEncoderInput,filter=1024, kernel_size=3, activation="relu", padding="same")
leftEncoderLayer2 <- layer_conv_1d(leftEncoderLayer1,filter=516, kernel_size=3, activation="relu", padding="same")
leftEncoderLayer3 <- layer_flatten(leftEncoderLayer2)
leftEncoderLayer4 <- layer_dense(leftEncoderLayer3, units = 128, activation=activation.1)


# Right encoder (st)
rightEncoderInput <-  layer_input(shape=c(3000))
rightEncoderLayer1 <- layer_dense(rightEncoderInput, units = 516, activation=activation.1)
rightEncoderLayer2 <- layer_dense(rightEncoderLayer1, units = 128, activation=activation.1)



#test.leftEncoderInput <-  layer_input(shape=c(3000, 5), dtype="float32")
#test.leftEncoderLayer1 <- layer_conv_1d(test.leftEncoderInput,filter=1024, kernel_size=3, activation="relu", padding="same")


#Encoder Merge

encoderMergeLayer = layer_dense(units=128, activation=activation.1)
leftMerge = encoderMergeLayer(leftEncoderLayer4)
rightMerge = encoderMergeLayer(rightEncoderLayer2)

#Different Merged Layer 

mergedLayer <-  layer_average(list(leftMerge, rightMerge))
leftMergedLayer <-  layer_average(list(leftMerge, leftMerge))
rightMergedLayer <-  layer_average(list(rightMerge, rightMerge))

z_mean <-  layer_dense(units = 16 )
z_log_sigma <-  layer_dense(units = 16)


# These three sets are used in differen models

sampling <- function(arg){
  z_mean <- arg[, 1:(16)]
  z_log_sigma <- arg[, (16 + 1):(2 * 16)]
  epsilon <- k_random_normal(shape = c(k_shape(z_mean)[[1]]), mean=0., stddev=epsilon_std)
  z_mean + k_exp(z_log_sigma)*epsilon
}



z_meanLeft <-  z_mean(leftMergedLayer)
z_log_sigmaLeft <-  z_log_sigma(leftMergedLayer)

z_meanRight <-  z_mean(rightMergedLayer)
z_log_sigmaRight <-  z_log_sigma(rightMergedLayer)

z_meanFull <-  z_mean(mergedLayer)
z_log_sigmaFull <-  z_log_sigma(mergedLayer)


zLeft <- 
  layer_concatenate(list(z_meanLeft, z_log_sigmaLeft)) %>% 
  layer_lambda(sampling)

zRight <- 
  layer_concatenate(list(z_meanRight, z_log_sigmaRight)) %>% 
  layer_lambda(sampling)
zFull <- 
  layer_concatenate(list(z_meanFull, z_log_sigmaFull)) %>% 
  layer_lambda(sampling)


leftEncoder = keras_model(leftEncoderInput, zLeft)
rightEncoder = keras_model(rightEncoderInput, zRight)
self.fullEncoder = keras_model(list(leftEncoderInput, rightEncoderInput), zFull)

# Left Dencoder (scData)
decoderInputs <-  layer_input(shape=(16))
decoderFirstLayer = layer_dense(decoderInputs, units=16, activation='relu')

leftdecoderLayer2 = layer_dense(decoderFirstLayer, units=c(5*3000),activation='relu')
leftdecoderLayer3 = layer_reshape(leftdecoderLayer2, target_shape = c(3000,5))
leftdecoderLayer4 = layer_conv_1d(leftdecoderLayer3,filter=516, kernel_size=3, activation="relu", padding="same")
leftdecoderLayer5 = layer_conv_1d(leftdecoderLayer4,filter=1024, kernel_size=3, activation="relu", padding="same")
leftdecoderLayerout = layer_conv_1d(leftdecoderLayer5,filter=5, kernel_size=3, activation="relu", padding="same")

# right Dencoder (stData)
rightdecoderLayer1 = layer_dense(decoderFirstLayer, units=256, activation=activation.1)
rightdecoderLayer2 = layer_dense(rightdecoderLayer1, units=516, activation=activation.1)
rightdecoderLayerout = layer_dense(rightdecoderLayer2, units=3000, activation='sigmoid')


#Decoder models

self.fullDecoder <-  keras_model(decoderInputs, list(leftdecoderLayerout, rightdecoderLayerout))
leftDeocder <-  keras_model(decoderInputs, leftdecoderLayerout)
rightDecoder <-  keras_model(decoderInputs, rightdecoderLayerout)



# Left to Right transition
outputs <-  self.fullDecoder(leftEncoder(leftEncoderInput))
self.leftToRightModel <-  keras_model(leftEncoderInput, outputs)
# leftToRightModel %>% summary()

# Right to Left transition
outputs <-  self.fullDecoder(rightEncoder(rightEncoderInput))
self.rightToLeftModel <-  keras_model(rightEncoderInput, outputs)
# rightToLeftModel %>% summary()

# Full Model
outputs = self.fullDecoder(self.fullEncoder(list(leftEncoderInput, rightEncoderInput)))

# Create the full model
vae_model = keras_model(list(leftEncoderInput, rightEncoderInput), outputs)

lowLearnAdam <- keras::optimizer_adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999, epsilon=NULL, decay=0.0, amsgrad=F)

vae_model %>% compile(optimizer=lowLearnAdam,loss=vae_loss)  # Compile
vae_model %>% summary()





# Freeze all shared layers
leftMerge$trainable = F
rightMerge$trainable = F

mergedLayer$trainable = F
leftMergedLayer$trainable = F
rightMergedLayer$trainable = F

z_meanLeft$trainable = F
z_log_sigmaLeft$trainable = F

z_meanRight$trainable = F
z_log_sigmaRight$trainable = F

z_meanFull$trainable = F
z_log_sigmaFull$trainable = F

zLeft$trainable = F
zRight$trainable = F
zFull$trainable = F

decoderFirstLayer$trainable = F

# Left VAE model which can't train middle
outputs <-  leftDeocder(leftEncoder(leftEncoderInput))
self.leftModel <-  keras_model(leftEncoderInput, outputs)
self.leftModel %>% compile(optimizer=lowLearnAdam, loss=left_vae_loss)


# Right VAE model which can't train middle
outputs <-  rightDecoder(rightEncoder(rightEncoderInput))
self.rightModel <-  keras_model(rightEncoderInput, outputs)
self.rightModel %>% compile(optimizer=lowLearnAdam, loss=right_vae_loss)


# Make shared layers trainable
leftMerge$trainable = T
rightMerge$trainable = T

mergedLayer$trainable = T
leftMergedLayer$trainable = T
rightMergedLayer$trainable = T

z_meanLeft$trainable = T
z_log_sigmaLeft$trainable = T

z_meanRight$trainable = T
z_log_sigmaRight$trainable = T

z_meanFull$trainable = T
z_log_sigmaFull$trainable = T

zLeft$trainable = T
zRight$trainable = T
zFull$trainable = T

decoderFirstLayer$trainable = T

# Make separate layers frozen

leftdecoderLayer2$trainable = F
leftdecoderLayer3$trainable = F
leftdecoderLayer4$trainable = F
leftdecoderLayer5$trainable = F
leftdecoderLayerout$trainable = F

# right Dencoder (stData)
rightdecoderLayer1$trainable = F
rightdecoderLayer2$trainable = F
rightdecoderLayerout$trainable = F





# Define center model
outputs = self.fullDecoder(self.fullEncoder(list(leftEncoderInput, rightEncoderInput)))
# Create the center model
self.centerModel = keras_model(list(leftEncoderInput, rightEncoderInput), outputs)
self.centerModel %>% compile(optimizer=lowLearnAdam, loss=vae_loss)  # Compile



# Training ---------------------------------------------------------------



train_model <- function(leftDomain, rightDomain, epoch=20){
  leftDomain_noisy <-  leftDomain
  rightDomain_noisy <-  rightDomain
  #callback <- keras::EarlyStopping(monitor='loss',patience=3, verbose=0, mode='auto')
  
  vae_model %>% fit(list(leftDomain, rightDomain),
                    list(leftDomain, rightDomain),
                     epochs=epoch,
                     batch_size=25,
                     shuffle=T,
                     verbose=1)
  
  for(i in 1:epoch){
    message(paste("On EPOCH: ", i))
    
    self.centerModel %>% fit(list(leftDomain, rightDomain),
                             list(leftDomain, rightDomain),
                             epochs=1,
                             batch_size=25,
                             shuffle=T)
    
    self.leftModel %>% fit(leftDomain, 
                           leftDomain,
                           epochs=1,
                           batch_size=25)
    
    self.rightModel %>% fit(rightDomain, 
                            rightDomain,
                            epochs=1,
                            batch_size=25)
  }
  
}

rightDomain <- array_reshape(data$st, c(1000,3000))
leftDomain <- array_reshape(data$sc, c(1000,3000,5))
for(i in 1:1000){leftDomain[i,,] <- as.numeric(leftDomain[i,,]+1)} #Add psydocount
for(i in 1:1000){rightDomain[i,] <- as.numeric(rightDomain[i,]+1)} #Add psydocount



train_model(leftDomain=leftDomain, rightDomain = rightDomain)



#Save training from all models


#vae_model %>% save_model_weights_tf("vae_model.ckpt")
#self.centerModel %>% save_model_weights_tf("self.centerModel.ckpt")
#self.leftModel %>% save_model_weights_tf("self.leftModel .ckpt")
#self.rightModel %>% save_model_weights_tf("self.rightModel.ckpt")
#saveRDS(rightDomain, "rightDomain.RDS")
#saveRDS(leftDomain, "leftDomain.RDS")




# Generator functions -----------------------------------------------------


# A|B -> A'|B'
pred <- vae_model %>% predict(list(leftDomain, rightDomain))

# Get integrated latent space


# Transfer A -> B'
# Predict st from sc 
rightDomain.pred <- self.leftToRightModel %>% predict(leftDomain)
rightDomain.pred %>% dim()

# Transfer  B -> A'
# Predict sc from st
leftDomain.pred <- self.rightToLeftModel %>% predict(rightDomain)
leftDomain.pred %>% dim()



# Validation both VAE -----------------------------------------------------
sc.pred.all <- pred[[1]]
st.pred.all <- pred[[2]]


sc.pred.all %>% dim()
sc.pred.all[1, 1:10, 1:5]

leftDomain[1, 1:10, 1:5]


sc.pred.all %>% dim()
st.pred.all[1, 1:10]
rightDomain[1, 1:10]



# Validation sc to st -----------------------------------------------------

st.org <- rightDomain
st.pred <- rightDomain.pred[[2]]

# Format array to count matrix
st.org <- map(.x=1:1000, .f=function(i){
  out <- st.org[i,] %>% as.data.frame()
  names(out) <- paste0("Cell_",i)
  rownames(out) <- var.genes
  return(out)
}) %>% do.call(cbind, .) %>% as.data.frame()
st.org %>% dim()
st.org[1:10, 1:5]


st.pred <- map(.x=1:1000, .f=function(i){
  out <- st.pred[i,] %>% as.data.frame()
  names(out) <- paste0("Cell_pred_",i)
  rownames(out) <- var.genes
  return(out)
}) %>% do.call(cbind, .) %>% as.data.frame()
st.pred %>% dim()

st.pred[1:10, 1:5]

plot(st.pred[,1] , st.org[,1], xlim = c(0,50), ylim = c(0,50), pch=19)







# Validation st to sc  -------------------------------------------------------------


sc.org <- leftDomain
sc.pred <- rightDomain.pred[[1]]

# Format array to count matrix
sc.org.df <- map(.x=1:1000, .f=function(i){
  out <- sc.org[i,,] %>% as.data.frame()
  names(out) <- paste0("Cell_",i,"_",1:5)
  rownames(out) <- var.genes
  return(out)
}) %>% do.call(cbind, .) %>% as.data.frame()
sc.org.df %>% dim()
sc.org.df[1:10, 1:5]


sc.pred <- map(.x=1:1000, .f=function(i){
  out <- sc.pred[i,,] %>% as.data.frame()
  names(out) <- paste0("Cell_pred_",i,"_",1:5)
  rownames(out) <- var.genes
  return(out)
}) %>% do.call(cbind, .) %>% as.data.frame()
sc.pred %>% dim()

sc.pred[1:10, 1:5]
#Transforme to counts
for(i in 1:ncol(sc.pred)){sc.pred[, i] <- round(sc.pred[, i])}
sc.pred[sc.pred == 1] <- 0

#Test 1. Scale the individual genes be a softmax

softmax <- function(x){exp(x)/sum(exp(x))}
x=sc.pred[1,] %>% as.numeric()
plot(x)
base <- c(min(x), max(x))
x %>% softmax() %>% scales::rescale(c(0,1)) %>% plot()


## Correlations
plot(sc.pred[,1:50] %>% rowMeans(), sc.org.df[,1:50] %>% rowMeans(), xlim = c(0,100), ylim = c(0,100))
plot(sc.pred[,1] , sc.org.df[,1], xlim = c(0,50), ylim = c(0,50), pch=19)
points(sc.pred[,2] , sc.org.df[,3], xlim = c(0,50), ylim = c(0,50), col="red",pch=19)
points(sc.pred[,3] , sc.org.df[,4], xlim = c(0,50), ylim = c(0,50), col="green",pch=19)
points(sc.pred[,4] , sc.org.df[,4], xlim = c(0,50), ylim = c(0,50), col="blue",pch=19)
points(sc.pred[,5] , sc.org.df[,4], xlim = c(0,50), ylim = c(0,50), col="orange",pch=19)

sc.pred.test <- sc.pred[, 1:50]
for(i in 1:nrow(sc.pred.test)){
  if(sc.pred.test[i, ] %>% sum()==0){sc.pred.test[i, ]=sc.pred.test[i, ]}
  else{sc.pred.test[i, ]<- scale(as.numeric(sc.pred.test[i, ]))}}
sc.pred.test <- scale(sc.pred.test)
pca <- irlba::prcomp_irlba(sc.pred.test)
pca$rotation %>% plot()


sc.org.seurat <- Seurat::CreateSeuratObject(sc.org.df,project="org" )
sc.pred.seurat <- Seurat::CreateSeuratObject(sc.pred, project="pred" )


sc.org.seurat <- 
  sc.org.seurat %>% 
  Seurat::FindVariableFeatures() %>%
  Seurat::ScaleData() %>% 
  Seurat::RunPCA() %>% 
  Seurat::RunUMAP(dims=1:10)

sc.org.seurat@meta.data$type <- "org"
Seurat::DimPlot(sc.org.seurat)

sc.pred.seurat <- 
  sc.pred.seurat %>% 
  Seurat::FindVariableFeatures() %>%
  Seurat::ScaleData() %>% 
  Seurat::RunPCA() %>% 
  Seurat::RunUMAP(dims=1:3)
Seurat::DimPlot(sc.pred.seurat)
sc.pred.seurat@meta.data$type <- "pred"


#Test2 Integrate data

sc.pred.seurat@meta.data$Set <- paste0("Set_", rep(1:5, 100))
merged.pred <- SeuratWrappers::RunFastMNN(Seurat::SplitObject(sc.pred.seurat, split.by = "Set" ), features = 100)
merged.pred <- merged.pred %>% Seurat::RunUMAP(dims=1:5, reduction="mnn")
Seurat::DimPlot(merged.pred, reduction = "mnn", pt.size = 0.1)





merged <- SeuratWrappers::RunFastMNN(list(sc.org.seurat, sc.pred.seurat), features = 100)

merged@meta.data

merged <- 
  merged %>% 
  Seurat::FindVariableFeatures() %>%
  Seurat::ScaleData() %>% 
  Seurat::RunPCA() %>% 
  Seurat::RunUMAP(dims=1:10, reduction="mnn")
Seurat::DimPlot(merged)







# Optimize model by optimite AE for sc data -------------------------------



# AE that is used 

leftDomain %>% dim()



# Left Encoder (scData)
leftEncoderInput <-  layer_input(shape=c(3000, 5), dtype="float32")
leftEncoderLayer3 <- layer_flatten(leftEncoderInput)
leftEncoderLayer4 <- layer_dense(leftEncoderLayer3, units = 256, activation="relu")
leftEncoderLayer5 <- layer_dense(leftEncoderLayer4, units = 128, activation="relu")
bt <- layer_dense(leftEncoderLayer5, units = 16)
decoderFirstLayer = layer_dense(bt, units=128, activation='relu')
leftdecoderLayer1 <- layer_dense(decoderFirstLayer, units = 256, activation="relu")
leftdecoderLayer4 = layer_dense(leftdecoderLayer1, units=c(5*3000),activation='softmax')
leftdecoderLayerout = layer_reshape(leftdecoderLayer4, target_shape = c(3000,5))





test.model.AE.autoencoder <- 
  keras_model(inputs = leftEncoderInput, outputs = leftdecoderLayerout) %>% 
  compile(optimizer="adam", loss="mse") 


leftDomain %>% dim()
test.model.AE.autoencoder %>% fit(leftDomain[1:200,,], epoch=20, batch_size=25)

pred.test <- test.model.AE.autoencoder %>% predict(leftDomain[1:5,,])

pred.test[1,1:10,]
leftDomain[1,1:10,]








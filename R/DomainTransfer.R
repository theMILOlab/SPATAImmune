


# Script Domain Transfer in multiomic data for integrative clustering and multiomic prediction
# Subject of multiomic data: Gene expression and metabolic data 


# Load Data ---------------------------------------------------------------
# Data are gene expression matrix and metabolic intensities

library(SPATA2)
exp <- readRDS("~/Desktop/SpatialTranscriptomics/Visium/Visium/ImageCoregistration/275_T/stSubset_296.RDS")
met <- readRDS("~/Desktop/SpatialTranscriptomics/Visium/Visium/ImageCoregistration/275_T/Metabolic_data_296.RDS")


#Intensity matrix

exp.mat <- SPATA2::getExpressionMatrix(exp)
dim(exp.mat)

met.mat <- SPATA2::getExpressionMatrix(met)
dim(met.mat)


# Prepare model -----------------------------------------------------------

shape <- c(ncol(met.mat))

if (tensorflow::tf$executing_eagerly())
  tensorflow::tf$compat$v1$disable_eager_execution()

library(keras)
K <- keras::backend()
epsilon_std <- 1.0
n.spots=1000
activation.1="tanh"
activation.2="sigmoid"
activation.3="relu"
latent_dim=16

# Loss Functions  ---------------------------------------------------------

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



# Encoders ----------------------------------------------------------------

# Left Encoder (exp)
library(keras)
leftEncoderInput <-  layer_input(shape=nrow(exp.mat), dtype="float32")
leftEncoderLayer1 <- layer_dense(leftEncoderInput, units = 516, activation=activation.1)
leftEncoderLayer2 <- layer_dense(leftEncoderLayer1, units = 128, activation=activation.1)
leftEncoderLayer3 <- layer_dense(leftEncoderLayer2, units = 64, activation=activation.1)


# Right encoder (met)
rightEncoderInput <-  layer_input(shape=nrow(met.mat), dtype="float32")
rightEncoderLayer1 <- layer_dense(rightEncoderInput, units = 516, activation=activation.1)
rightEncoderLayer2 <- layer_dense(rightEncoderLayer1, units = 128, activation=activation.1)
rightEncoderLayer3 <- layer_dense(rightEncoderLayer2, units = 64, activation=activation.1)


#Encoder Merge
encoderMergeLayer = layer_dense(units=64, activation=activation.1)
leftMerge = encoderMergeLayer(leftEncoderLayer3)
rightMerge = encoderMergeLayer(rightEncoderLayer3)



# Latent Space ------------------------------------------------------------

#Different Merged Layer 
mergedLayer <-  layer_average(list(leftMerge, rightMerge))
leftMergedLayer <-  layer_average(list(leftMerge, leftMerge))
rightMergedLayer <-  layer_average(list(rightMerge, rightMerge))

z_mean <-  layer_dense(units = latent_dim )
z_log_sigma <-  layer_dense(units = latent_dim)


# These three sets are used in differen models
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
  layer_concatenate(list(z_meanFull, z_log_sigmaFull), name = "LS") %>% 
  layer_lambda(sampling)


leftEncoder = keras_model(leftEncoderInput, zLeft)
rightEncoder = keras_model(rightEncoderInput, zRight)
self.fullEncoder = keras_model(list(leftEncoderInput, rightEncoderInput), zFull)



# Decoder -----------------------------------------------------------------

# Left Dencoder (exp)
decoderInputs <-  layer_input(shape=(16))
decoderFirstLayer = layer_dense(decoderInputs, units=16, activation='relu')

leftdecoderLayer1 = layer_dense(decoderFirstLayer, units=64,activation='relu')
leftdecoderLayer2 = layer_dense(leftdecoderLayer1, units=128,activation='relu')
leftdecoderLayer3 = layer_dense(leftdecoderLayer2, units=516,activation='relu')
leftdecoderLayerout = layer_dense(leftdecoderLayer3, units=nrow(exp.mat),activation='relu')


# right Dencoder (met)
rightdecoderLayer1 = layer_dense(decoderFirstLayer, units=64,activation='relu')
rightdecoderLayer2 = layer_dense(rightdecoderLayer1, units=128,activation='relu')
rightdecoderLayer3 = layer_dense(rightdecoderLayer2, units=516,activation='relu')
rightdecoderLayerout = layer_dense(rightdecoderLayer3, units=nrow(met.mat),activation='relu')

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






# Freez Layer shared ------------------------------------------------------



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



# self Models -------------------------------------------------------------


# Left VAE model which can't train middle
outputs <-  leftDeocder(leftEncoder(leftEncoderInput))
self.leftModel <-  keras_model(leftEncoderInput, outputs)
self.leftModel %>% compile(optimizer=lowLearnAdam, loss=left_vae_loss)


# Right VAE model which can't train middle
outputs <-  rightDecoder(rightEncoder(rightEncoderInput))
self.rightModel <-  keras_model(rightEncoderInput, outputs)
self.rightModel %>% compile(optimizer=lowLearnAdam, loss=right_vae_loss)


# Trainable ---------------------------------------------------------------


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

# Seperate Frozen ---------------------------------------------------------
# Make separate layers frozen

leftdecoderLayer2$trainable = F
leftdecoderLayer3$trainable = F
leftdecoderLayer4$trainable = F
leftdecoderLayer5$trainable = F
leftdecoderLayerout$trainable = F
rightdecoderLayer1$trainable = F
rightdecoderLayer2$trainable = F
rightdecoderLayerout$trainable = F

# Compile Model -------------------------------------------------------------------

# Define center model
outputs = self.fullDecoder(self.fullEncoder(list(leftEncoderInput, rightEncoderInput)))
# Create the center model
self.centerModel = keras_model(list(leftEncoderInput, rightEncoderInput), outputs)
self.centerModel %>% compile(optimizer=lowLearnAdam, loss=vae_loss)  # Compile











# Training ----------------------------------------------------------------
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

leftDomain <- array_reshape(exp.mat, c(dim(exp.mat)[2:1]))
rightDomain <- array_reshape(met.mat, c(dim(met.mat)[2:1]))

class(rightDomain)
dim(rightDomain)
dim(leftDomain)

train_model(leftDomain, rightDomain, epoch=10)




# Generator functions -----------------------------------------------------


# A|B -> A'|B'
pred <- vae_model %>% predict(list(leftDomain, rightDomain))

# Transfer A -> B'
rightDomain.pred <- self.leftToRightModel %>% predict(leftDomain)
# Transfer  B -> A'
leftDomain.pred <- self.rightToLeftModel %>% predict(rightDomain)

get.latent.space <- function(leftDomain, rightDomain){
  keras_model(inputs = vae_model$input, 
              outputs = get_layer(vae_model, "model_89")$output) %>% predict(list(leftDomain, rightDomain))
  
} 
LS <- get.latent.space(leftDomain, rightDomain)


#Latent space
umap <- uwot::umap(LS,n_neighbors = 5, learning_rate = 0.5, init = "pca", n_epochs = 20)
umap <- umap %>% as.data.frame()
names(umap) <- c("UMAP1", "UMAP2")
umap %>% 
  #mutate(class = as.factor(mnist$test$y)) %>%
  ggplot(aes(x = UMAP1, y = UMAP2)) + 
  geom_point()+
  theme_classic()

leftDomain %>% dim()
rightDomain.pred[[1]] %>% dim()


rightDomain.pred[[2]] %>% dim()
rightDomain %>% dim()


coords <- SPATA2::getCoordsDf(exp) %>% as.data.frame()


#Check Metabolites

z=1340

plot(rightDomain.pred[[2]][,z], rightDomain[,z])

p1 <- ggplot(coords, aes(x = x, y = y, color=rightDomain.pred[[2]][,z])) + 
  geom_point(size=2)+
  theme_classic()+scale_color_viridis_c()
p2 <- ggplot(coords, aes(x = x, y = y, color=rightDomain[,z])) + 
  geom_point(size=2)+
  theme_classic()+scale_color_viridis_c()

library(patchwork)
p1+p2



pred.exp <- rightDomain.pred[[1]]
rownames(pred.exp) <- rownames(exp.mat %>% t())
colnames(pred.exp) <- colnames(exp.mat %>% t())
pred.exp <- t(pred.exp)

pred.exp %>% rowSums() 
gene="GFAP"


plot(pred.exp[gene, ], exp.mat[gene, ])

p1 <- ggplot(coords, aes(x = x, y = y, color=pred.exp[gene, ])) + 
  geom_point(size=3)+
  theme_classic()+scale_color_viridis_c()
p2 <- ggplot(coords, aes(x = x, y = y, color=exp.mat[gene, ])) + 
  geom_point(size=3)+
  theme_classic()+scale_color_viridis_c()

library(patchwork)
p1+p2


















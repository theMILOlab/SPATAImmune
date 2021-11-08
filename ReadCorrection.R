



IgBlast <- read.table("middle.tsv", sep="\t", header=T)


#get training data set 
library(tidyverse)
IgBlast.f %>% head()


library(Biostrings)
IgBlast.f <- 
  IgBlast %>% 
  filter(cdr3_aa != "")

IgBlast.f$cdr3= Biostrings::translate(IgBlast.f$cdr3_aa %>% DNAStringSet()) %>% as.data.frame() %>% pull(x)
dim(IgBlast.f)

# Check locus 
# alpha with V and J segment, ß with v d and J segment 


vdj <- IgBlast.f %>% select(locus, d_call,j_call,sequence_alignment)
names(vdj) <-c("c", "v", "d", "j") 
vdj$number=1:nrow(vdj)
neg.select <- vdj %>% filter(c=="TRB" & d=="") %>% pull(number)
vdj <- vdj[!vdj$number %in% neg.select, ]
IgBlast.f <- IgBlast.f[vdj$number, ]


vdj_TRB <- vdj %>% filter(c=="TRB") %>% mutate(reAr=paste0(v,"_",d,"_",j))
vdj_TRB %>% group_by(reAr) %>% count() %>% arrange(desc(n))








#start with only vj in TRA

IgBlast.use <- IgBlast.f
IgBlast.f <- IgBlast.use %>% filter(locus == "TRA")
IgBlast.f %>% dim()

# Network to predict the VDJ rearrangment 
#1. Create overlaping kmers from the sequences with a total number of mx/4 and kmer length 4

CreateKmers <- function(seq, kmer.length=10, kmer.nr=NULL, overlap=2){
  
  length.seq <- stringi::stri_length(seq)
  max.seq.length <- max(length.seq, na.rm = T)
  min.seq.length <- min(length.seq, na.rm = T)
  
  if(is.null(kmer.nr)){
    kmer.nr=min.seq.length/4
  }else{kmer.nr=kmer.nr}
  message(paste0(Sys.time(), " --- Sequence range between: ",min.seq.length, " to ", max.seq.length, " bp ---" ))
  
  kmers <- map(.x=seq, .f=function(s){
    
    seq.x <- unlist(strsplit(s, "") )
    nr.kmers <- 1+(length(seq.x)/(kmer.length-overlap)) %>% round()
    v=seq(1,(length(seq.x)-kmer.length), by=8)
    x <- map(.x=v, .f=function(i){seq.x[i:(i+kmer.length)] %>% paste0(collapse="")}) %>% unlist() %>% paste(collapse = " ")
    
    
  }) %>% unlist()

}
  
kmers <- CreateKmers(seq=IgBlast.f$sequence[1:100000], kmer.length=5)



kmers <- iconv(kmers)
head(kmers, 2)

library(reticulate)
library(purrr)
skipgrams_generator <- function(kmers, tokenizer, window_size, negative_samples) {
  gen <- texts_to_sequences_generator(tokenizer, sample(kmers))
  function() {
    skip <- generator_next(gen) %>%
      skipgrams(
        vocabulary_size = tokenizer$num_words, 
        window_size = window_size, 
        negative_samples = 1
      )
    x <- transpose(skip$couples) %>% map(. %>% unlist %>% as.matrix(ncol = 1))
    y <- skip$labels %>% as.matrix(ncol = 1)
    list(x, y)
  }
}

# Create Tokenizer
tokenizer <- 
  text_tokenizer(200000) %>% 
  fit_text_tokenizer(kmers)
tokenizer <- 
  text_tokenizer(tokenizer$word_index %>% length()) %>% fit_text_tokenizer(kmers)

#Parameter
embedding_size <- 200 
skip_window <- 5
num_sampled <- 1


word2vec <- function(){

library(keras)

input_target <- layer_input(shape = 1)
input_context <- layer_input(shape = 1)
  
embedding <- layer_embedding(
    input_dim = tokenizer$num_words+1, 
    output_dim = embedding_size, 
    input_length = 1, 
    name = "embedding")
  
target_vector <- 
    input_target %>% 
    embedding() %>% 
    layer_flatten(name = "Test")


context_vector <- 
    input_context %>%
    embedding() %>%
    layer_flatten()

  
dot_product <- layer_dot(list(target_vector, context_vector), axes = 1)
output <- layer_dense(dot_product, units = 1, activation = "sigmoid")
  
model <- keras_model(list(input_target, input_context), output)
model %>% compile(loss = "binary_crossentropy", optimizer = "adam")

return(model)

} 

model <- word2vec()

model %>% 
  fit_generator(
  skipgrams_generator(kmers, tokenizer, skip_window, negative_samples=1), 
  steps_per_epoch = 10, 
  epochs = 40)



# Get vector presentation


embedding_matrix <- get_weights(model)[[1]]
embedding_matrix %>% dim()


#plot(umap::umap(embedding_matrix)$layout)

words <- data_frame(
  word = names(tokenizer$word_index), 
  id = as.integer(unlist(tokenizer$word_index))
)

words <- words %>%
  filter(id <= tokenizer$num_words) %>%
  arrange(id)

row.names(embedding_matrix) <- c("UNK", words$word)




library(text2vec)
find_similar_words <- function(word, embedding_matrix, n = 5) {
  similarities <- embedding_matrix[word, , drop = FALSE] %>%
    sim2(embedding_matrix, y = ., method = "cosine")
  
  similarities[,1] %>% sort(decreasing = TRUE) %>% head(n)
}
find_similar_words("ggtatcaacgc", embedding_matrix)



#Represent sequences by the kmers embedded matrix

#Test LSTM for sequence embedding

# set some parameters for our model
# the number of kmers per sequence
max_len <- map(.x=kmers[1:1000], .f=function(k){
  length=str_split(k, pattern=" ") %>% unlist %>% length()
}) %>% unlist() %>% max()  


# Use 0 padding to avarage shape
list2ary = function(input.list){  #input a list of lists
  rows.cols <- dim(input.list[[1]])
  sheets <- length(input.list)
  output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
  colnames(output.ary) <- colnames(input.list[[1]])
  row.names(output.ary) <- row.names(input.list[[1]])
  return(output.ary)    # output as a 3-D array
}
embedding_seq <- map(.x=kmers[1:1000], .f=function(k){
  mat <- embedding_matrix[str_split(k, pattern=" ") %>% unlist %>%  tolower(), ] 
  mat <- rbind(mat, matrix(0, max_len-nrow(mat), ncol(mat)))
  })
array <- list2ary(embedding_seq) %>% aperm(., c(3,1,2))
dim(array)


# Parameter ---------------------------------------------------------------


if (tensorflow::tf$executing_eagerly())
  tensorflow::tf$compat$v1$disable_eager_execution()
library(keras)
K <- keras::backend()
timesteps <- max_len



# LSTM Autoencoder --------------------------------------------------------


model.LSTM.encoder <- 
  keras_model_sequential() %>% 
  layer_masking(input_shape=c(timesteps, 200), mask_value=0) %>% 
  layer_lstm(units=254, input_shape=c(timesteps, 200), activation="relu", return_sequences = T) %>% 
  layer_lstm(units=128, activation="relu") %>% 
  layer_dense(units = 128, name ="latent_space") %>% 
  layer_repeat_vector(timesteps) %>% 
  layer_lstm(units=timesteps ,return_sequences=T, activation="relu") %>% 
  time_distributed(layer = layer_dense(units = 200))


model.LSTM.encoder  %>% compile(loss = 'mae', optimizer = 'adam')

model.LSTM.encoder %>% fit(array, array, epochs=10)


#get latent space
get.latent.space.2D.LSTM.Autoencoder <- function(model,array){
  keras_model(inputs = model$input, 
              outputs = get_layer(model, "latent_space")$output) %>% 
    predict(., array)
  
} 

ls <- get.latent.space.2D.LSTM.Autoencoder(model=model.LSTM.encoder, array)
ls %>% dim()

plot(umap::umap(ls)$layout)



# Test Transformers -------------------------------------------------------






# Prediction --------------------------------------------------------------


# Now sequences are embedded into a fixed n*128 dimensional space
# Now predict the alignments and CD3 Sequence


# transforme the CDR3 regions into kmers followed by padding and transformation using the word2vec library
#













# Predict the Vdj arrangment based on the latent space 

J.segment <- IgBlast.f$j_family[1:1000]

dim(ls)

Jx_train <- ls[1:500, ]
Jy_train <- J.segment %>% as.factor() %>% as.numeric() %>% to_categorical()
colnames(Jy_train[,2:3]) <- c(unique(J.segment))
Jy_train <- Jy_train[,2:3]


#Model 

model.J <- 
  keras::keras_model_sequential() %>% 
  keras::layer_dense(20, input_shape = 128) %>% 
  keras::layer_dropout(0.5) %>% 
  keras::layer_dense(10,activation='relu') %>% 
  keras::layer_dropout(0.2) %>% 
  keras::layer_dense(2, activation='softmax')

model.J  %>% compile(loss = 'categorical_crossentropy', optimizer = 'adam')
model.J %>% fit(ls[1:1000, ], Jy_train[1:1000, ], epoch=20)


#Predict

pred <- predict(model.J, ls[101:200, ]) %>% as.data.frame()
pred$J.segment <- IgBlast.f$j_family[101:200]



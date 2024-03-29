% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/family.run.R
\name{runSPTCRdeep}
\alias{runSPTCRdeep}
\title{runSPTCRdeep -Deprecated -}
\usage{
runSPTCRdeep(
  object,
  constant = NULL,
  run.AA2vector = T,
  run.ANN = T,
  kmer.length = 3,
  embedding_size = 10,
  skip_window = 5,
  epochs.w2v = 80,
  LSTM = T,
  epochs.LSTM = 40,
  VAE = F,
  num_sampled = 1,
  epochs = 40,
  LSTMBatch = 1000,
  temp.save = T,
  temp.folder = getwd()
)
}
\value{
SPATA object
}
\description{
This function import the VDJ data from the SPTCR-Pipline or IgBlast
}
\author{
Dieter Henrik Heiland
}

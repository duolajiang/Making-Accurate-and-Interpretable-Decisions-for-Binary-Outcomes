library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))
library(BART)
source('simulation.R')
source('main_function.R')

# scenario <- 'H'
# true.names.Z1 <- ''
# true.value.Z1 <- ''
# true.names.Z0 <- c("Ca","Cb")
# true.value.Z0 <- c(0.5,0.5)
# n_pop <- 10000
# n_sample <- 1000
# n_sim <- 100
# complexity <- 10
# FUN <- estimation_bart_tree
# lambda <- log(3)
# file_name='BART_TREE'
# ntree <- 200
# 
# SimulationMain_variance_bias_PEHE (scenario = 'H',
#                                    true.names.Z1 = true.names.Z1,
#                                    true.value.Z1 = true.value.Z1,
#                                    true.names.Z0 = true.names.Z0 ,
#                                    true.value.Z0 = true.value.Z0,
#                                    n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
#                                    complexity=complexity,
#                                    FUN=FUN,lambda=lambda,
#                                    file_name=file_name,ntree=ntree)


scenario <- 'H'
true.names.Z1 <- ''
true.value.Z1 <- ''
true.names.Z0 <- c("Ca","Cb")
true.value.Z0 <- c(0.5,0.5)
n_pop <- 10000
n_sample <- 1000
n_sim <- 100
FUN <- estimation_bart_lreg
lambda <- log(3)
eta <- 0.01
alpha <- 1
step <- 0.05
epoch <- 1000
batch <- 1
ntree <- 200
L2 <- FALSE
ct00 <- 1
ct01 <- 1
ct10 <- 0
ct11 <- 0
cc00 <- 1
cc01 <- 0
cc10 <- 1
cc11 <- 0
rho <- 0
file_name=glue::glue('BART_LReg_',ct10)

SimulationMain_LReg (scenario = 'H',
                     true.names.Z1 = true.names.Z1,
                     true.value.Z1 = true.value.Z1,
                     true.names.Z0 = true.names.Z0 ,
                     true.value.Z0 = true.value.Z0,
                     n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
                     alpha=alpha,step=step,epoch=epoch,eta=eta,batch=batch,
                     FUN=FUN,
                     file_name=file_name,ntree=ntree,L2=L2,
                     ct00=ct00, ct01=ct01,ct10=ct10,ct11=ct11,
                     cc00=cc00, cc01=cc01,cc10=cc10,cc11=cc11,
                     rho)

ct00 <- 1.05
ct01 <- 1.05
ct10 <- 0.05
ct11 <- 0.05
file_name=glue::glue('BART_LReg_',ct10)
SimulationMain_LReg (scenario = 'H',
                     true.names.Z1 = true.names.Z1,
                     true.value.Z1 = true.value.Z1,
                     true.names.Z0 = true.names.Z0 ,
                     true.value.Z0 = true.value.Z0,
                     n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
                     alpha=alpha,step=step,epoch=epoch,eta=eta,batch=batch,
                     FUN=FUN,
                     file_name=file_name,ntree=ntree,L2=L2,
                     ct00=ct00, ct01=ct01,ct10=ct10,ct11=ct11,
                     cc00=cc00, cc01=cc01,cc10=cc10,cc11=cc11,
                     rho)

ct00 <- 1.1
ct01 <- 1.1
ct10 <- 0.1
ct11 <- 0.1
file_name=glue::glue('BART_LReg_',ct10)
SimulationMain_LReg (scenario = 'H',
                     true.names.Z1 = true.names.Z1,
                     true.value.Z1 = true.value.Z1,
                     true.names.Z0 = true.names.Z0 ,
                     true.value.Z0 = true.value.Z0,
                     n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
                     alpha=alpha,step=step,epoch=epoch,eta=eta,batch=batch,
                     FUN=FUN,
                     file_name=file_name,ntree=ntree,L2=L2,
                     ct00=ct00, ct01=ct01,ct10=ct10,ct11=ct11,
                     cc00=cc00, cc01=cc01,cc10=cc10,cc11=cc11,
                     rho)

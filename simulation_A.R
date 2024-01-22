library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))
library(BART)
source('simulation.R')
source('main_function.R')

# scenario <- 'A'
# true.names.Z1 <- c("C1", "B1", "a3:A1")
# true.value.Z1 <- c(0.5, 2, 2)
# true.names.Z0 <- c("C1")
# true.value.Z0 <- c(0.5)
# n_pop <- 10000
# n_sample <- 1000
# n_sim <- 100
# complexity <- 10
# FUN <- estimation_bart_tree
# lambda <- log(3)
# file_name='BART_TREE'
# ntree <- 200
# 
# #simulate m1: BART, m2:BART-TREE, m3:TREE
# SimulationMain_variance_bias_PEHE (scenario = 'A',
#                                    true.names.Z1 = true.names.Z1,
#                                    true.value.Z1 = true.value.Z1,
#                                    true.names.Z0 = true.names.Z0 ,
#                                    true.value.Z0 = true.value.Z0,
#                                    n_pop=n_pop,n_sample=n_sample, n_sim=n_sim,
#                                    complexity=complexity,
#                                    FUN=FUN,lambda=lambda,
#                                    file_name=file_name,ntree=ntree)



scenario <- 'A'
true.names.Z1 <- c("C1", "B1", "a3:A1")
true.value.Z1 <- c(0.5, 2, 2)
true.names.Z0 <- c("C1")
true.value.Z0 <- c(0.5)
n_pop <- 10000
n_sample <- 1000
n_sim <- 100
FUN <- estimation_bart_lreg
lambda <- log(3)
eta <- 0.01
alpha <- 1
step <- 0.05
epoch <- 1 # number of loops for updating weights using the training dataset
batch <- 1 #actually lead to n_sample/batch iterations in each epoch, when batch=1, it is SGD. Combined with epoch, the total number of iterations for updating weights is epoch*n_sample/batch
ntree <- 200
L2 <- FALSE
cc00 <- 1
cc01 <- 0
cc10 <- 1
cc11 <- 0
rho <- 0

for(loss_treatment in seq(0.2,0.4,0.05)){
  ct00 <- 1+loss_treatment
  ct01 <- 1+loss_treatment
  ct10 <- loss_treatment
  ct11 <- loss_treatment
  file_name=glue::glue('BART_LReg_',ct10)
  SimulationMain_LReg (scenario = 'A',
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
  
}

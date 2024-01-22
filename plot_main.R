library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))
source("plot_function.R")
library(ggplot2)
library(ggpubr)

# =======================================================
# plot loss outcome by model complexity; generate plots for each complexity parameter. 
# in each model complexity plot, 8 plots for each scenario, and x-axis is threhold, y is outcome/loss in each subplot. 
# =======================================================
threhold.max <- 0.4
step <- 0.05
filename <- 'results/BART_TREE10000_'
for (depth in seq(10)){
  load(paste(filename,'A.RData',sep='')); p.A <- plot_BART_TREE_loss_outcome(threhold.max, step, depth, outcome)
  load(paste(filename,'B.RData',sep='')); p.B <- plot_BART_TREE_loss_outcome(threhold.max, step, depth, outcome)
  load(paste(filename,'C.RData',sep='')); p.C <- plot_BART_TREE_loss_outcome(threhold.max, step, depth, outcome)
  load(paste(filename,'D.RData',sep='')); p.D <- plot_BART_TREE_loss_outcome(threhold.max, step, depth, outcome)
  load(paste(filename,'E.RData',sep='')); p.E <- plot_BART_TREE_loss_outcome(threhold.max, step, depth, outcome)
  load(paste(filename,'F.RData',sep='')); p.F <- plot_BART_TREE_loss_outcome(threhold.max, step, depth, outcome)
  load(paste(filename,'G.RData',sep='')); p.G <- plot_BART_TREE_loss_outcome(threhold.max, step, depth, outcome)
  load(paste(filename,'H.RData',sep='')); p.H <- plot_BART_TREE_loss_outcome(threhold.max, step, depth, outcome)
  
  pdf(file=glue::glue("plots/simulation_loss_outcome_TREE_",depth,".pdf"),width = 16,height = 15,onefile=FALSE) 
  print(ggarrange(
    p.A[[1]], p.B[[1]], p.C[[1]], p.D[[1]],
    p.E[[1]], p.F[[1]], p.G[[1]], p.H[[1]],
    p.A[[2]], p.B[[2]], p.C[[2]], p.D[[2]],
    p.E[[2]], p.F[[2]], p.G[[2]], p.H[[2]],
    labels=c(
      'A','B','C','D',
      'E','F','G','H',
      'A','B','C','D',
      'E','F','G','H'),
    nrow = 4,ncol=4,common.legend = TRUE,legend = 'bottom'))
  dev.off() 
}

# # =======================================================
# # plot precision and recall when threhold=0
# # generate a plot with 8 subplots;
# # in each subplot, x-axis is model complexity (i.e., depth), y is precision/recall
# # =======================================================
#complexity <- 10
#filename <- 'results/BART_TREE10000_'
#threhold <- 0
# load(paste(filename,'A.RData',sep='')); p.A <- plot_accuracy_metrics(outcome,complexity,threhold)
# load(paste(filename,'B.RData',sep='')); p.B <- plot_accuracy_metrics(outcome,complexity,threhold)
# load(paste(filename,'C.RData',sep='')); p.C <- plot_accuracy_metrics(outcome,complexity,threhold)
# load(paste(filename,'D.RData',sep='')); p.D <- plot_accuracy_metrics(outcome,complexity,threhold)
# load(paste(filename,'E.RData',sep='')); p.E <- plot_accuracy_metrics(outcome,complexity,threhold)
#load(paste(filename,'F.RData',sep='')); p.F <- plot_accuracy_metrics(outcome,complexity,threhold)
# load(paste(filename,'G.RData',sep='')); p.G <- plot_accuracy_metrics(outcome,complexity,threhold)
# load(paste(filename,'H.RData',sep='')); p.H <- plot_accuracy_metrics(outcome,complexity,threhold)
#   
# pdf(file=glue::glue("plots/simulation_precision_recall_threhold_",threhold,".pdf"),width = 11,height = 11,onefile=FALSE) 
# print(ggarrange(
#   p.A[[3]],p.B[[3]],p.C[[3]],p.D[[3]], # precision
#   p.E[[3]],p.F[[3]],p.G[[3]],p.H[[3]],
#   p.A[[4]],p.B[[4]],p.C[[4]],p.D[[4]], # recall
#   p.E[[4]],p.F[[4]],p.G[[4]],p.H[[4]],
#   labels=c(
#     'A','B','C','D','E','F','G','H',
#     'A','B','C','D','E','F','G','H'),
#   nrow = 4,ncol=4,common.legend = TRUE,legend = 'bottom'))
# dev.off() 
# 
# pdf(file=glue::glue("plots/simulation_tree_acc_",threhold,".pdf"),width = 11,height = 11,onefile=FALSE) 
# print(ggarrange(
#   p.A[[1]],p.B[[1]],p.C[[1]],p.D[[1]], # precision
#   p.E[[1]],p.F[[1]],p.G[[1]],p.H[[1]],
#   labels=c('A','B','C','D','E','F','G','H'),
#   nrow = 2,ncol=4,common.legend = TRUE,legend = 'bottom'))
# dev.off() 
# 
#   
# =======================================================
# plot accuracy, precision, recall when threhold 0-0.4
# =======================================================
threhold.max <- 0.4
step <- 0.05
filename <- 'results/BART_TREE10000_'
depth <- 2
load(paste(filename, 'A.RData', sep = ''))
p.A <- plot_BART_TREE_acc_precision_recall(threhold.max, step, depth, outcome)
load(paste(filename, 'B.RData', sep = ''))
p.B <- plot_BART_TREE_acc_precision_recall(threhold.max, step, depth, outcome)
load(paste(filename, 'C.RData', sep = ''))
p.C <- plot_BART_TREE_acc_precision_recall(threhold.max, step, depth, outcome)
load(paste(filename, 'D.RData', sep = ''))
p.D <- plot_BART_TREE_acc_precision_recall(threhold.max, step, depth, outcome)
load(paste(filename, 'E.RData', sep = ''))
p.E <- plot_BART_TREE_acc_precision_recall(threhold.max, step, depth, outcome)
load(paste(filename, 'F.RData', sep = ''))
p.F <- plot_BART_TREE_acc_precision_recall(threhold.max, step, depth, outcome)
load(paste(filename, 'G.RData', sep = ''))
p.G <- plot_BART_TREE_acc_precision_recall(threhold.max, step, depth, outcome)
load(paste(filename, 'H.RData', sep = ''))
p.H <- plot_BART_TREE_acc_precision_recall(threhold.max, step, depth, outcome)

pdf(
  file = glue::glue(
    "plots/simulation_BART_TREE_acc_precision_recall_depth_",
    depth,
    ".pdf"
  ),
  width = 15,
  height = 20,
  onefile = FALSE
)
print(
  ggarrange(
    # accuracy
    p.A[[1]],
    p.B[[1]],
    p.C[[1]],
    p.D[[1]],
    p.E[[1]],
    p.F[[1]],
    p.G[[1]],
    p.H[[1]],
    # precision
    p.A[[2]],
    p.B[[2]],
    p.C[[2]],
    p.D[[2]],
    p.E[[2]],
    p.F[[2]],
    p.G[[2]],
    p.H[[2]],
    # recall
    p.A[[3]],
    p.B[[3]],
    p.C[[3]],
    p.D[[3]],
    p.E[[3]],
    p.F[[3]],
    p.G[[3]],
    p.H[[3]], 
    labels = c(
      'A',
      'B',
      'C',
      'D',
      'E',
      'F',
      'G',
      'H',
      'A',
      'B',
      'C',
      'D',
      'E',
      'F',
      'G',
      'H',
      'A',
      'B',
      'C',
      'D',
      'E',
      'F',
      'G',
      'H'
    ),
    nrow = 6,
    ncol = 4,
    common.legend = TRUE,
    legend = 'bottom'
  )
)
dev.off() 


# =======================================================
# 
# =======================================================
# threhold.max <- 0.4
# step <- 0.05
# filename <- 'results/BART_TREE10000_'
# depth <- 2
# ct00 <- 1
# ct01 <- 1
# ct10 <- 0
# ct11 <- 0
# cc00 <- 1
# cc01 <- 0
# cc10 <- 1
# cc11 <- 0
# rho <- 0
# 
# load(paste(filename, 'A.RData', sep = ''))
# p.A <- plot_BART_TREE_acc_precision_recall_m1_from_BARTLReg_simulation(threhold.max, step, depth, outcome,filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "A",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
# load(paste(filename, 'B.RData', sep = ''))
# p.B <- plot_BART_TREE_acc_precision_recall_m1_from_BARTLReg_simulation(threhold.max, step, depth, outcome,filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "B",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
# load(paste(filename, 'C.RData', sep = ''))
# p.C <- plot_BART_TREE_acc_precision_recall_m1_from_BARTLReg_simulation(threhold.max, step, depth, outcome,filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "C",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
# load(paste(filename, 'D.RData', sep = ''))
# p.D <- plot_BART_TREE_acc_precision_recall_m1_from_BARTLReg_simulation(threhold.max, step, depth, outcome,filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "D",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
# load(paste(filename, 'E.RData', sep = ''))
# p.E <- plot_BART_TREE_acc_precision_recall_m1_from_BARTLReg_simulation(threhold.max, step, depth, outcome,filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "E",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
# load(paste(filename, 'F.RData', sep = ''))
# p.F <- plot_BART_TREE_acc_precision_recall_m1_from_BARTLReg_simulation(threhold.max, step, depth, outcome,filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "F",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
# load(paste(filename, 'G.RData', sep = ''))
# p.G <- plot_BART_TREE_acc_precision_recall_m1_from_BARTLReg_simulation(threhold.max, step, depth, outcome,filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "G",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
# load(paste(filename, 'H.RData', sep = ''))
# p.H <- plot_BART_TREE_acc_precision_recall_m1_from_BARTLReg_simulation(threhold.max, step, depth, outcome,filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "H",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
# 
# pdf(
#   file = glue::glue(
#     "plots/simulation_BART_TREE_acc_precision_recall_depth_",
#     depth,
#     ".pdf"
#   ),
#   width = 15,
#   height = 20,
#   onefile = FALSE
# )
# print(
#   ggarrange(
#     # accuracy
#     p.A[[1]],
#     p.B[[1]],
#     p.C[[1]],
#     p.D[[1]],
#     p.E[[1]],
#     p.F[[1]],
#     p.G[[1]],
#     p.H[[1]],
#     # precision
#     p.A[[2]],
#     p.B[[2]],
#     p.C[[2]],
#     p.D[[2]],
#     p.E[[2]],
#     p.F[[2]],
#     p.G[[2]],
#     p.H[[2]],
#     # recall
#     p.A[[3]],
#     p.B[[3]],
#     p.C[[3]],
#     p.D[[3]],
#     p.E[[3]],
#     p.F[[3]],
#     p.G[[3]],
#     p.H[[3]], 
#     labels = c(
#       'A',
#       'B',
#       'C',
#       'D',
#       'E',
#       'F',
#       'G',
#       'H',
#       'A',
#       'B',
#       'C',
#       'D',
#       'E',
#       'F',
#       'G',
#       'H',
#       'A',
#       'B',
#       'C',
#       'D',
#       'E',
#       'F',
#       'G',
#       'H'
#     ),
#     nrow = 6,
#     ncol = 4,
#     common.legend = TRUE,
#     legend = 'bottom'
#   )
# )
# dev.off() 



# =======================================================  
# =======================================================
# BART-logistic_regression, 
# plot average loss and outcome for threhold from 0~0.4
# =======================================================
cc00 <- 1
cc01 <- 0
cc10 <- 1
cc11 <- 0
rho <- 0
threhold.max <- 0.4
step <- 0.05

p.A <- plot_BART_LReg_loss_outcome(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "A",cc00,cc01,cc10,cc11,rho)
p.B <- plot_BART_LReg_loss_outcome(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "B",cc00,cc01,cc10,cc11,rho)
p.C <- plot_BART_LReg_loss_outcome(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "C",cc00,cc01,cc10,cc11,rho)
p.D <- plot_BART_LReg_loss_outcome(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "D",cc00,cc01,cc10,cc11,rho)
p.E <- plot_BART_LReg_loss_outcome(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "E",cc00,cc01,cc10,cc11,rho)
p.F <- plot_BART_LReg_loss_outcome(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "F",cc00,cc01,cc10,cc11,rho)
p.G <- plot_BART_LReg_loss_outcome(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "G",cc00,cc01,cc10,cc11,rho)
p.H <- plot_BART_LReg_loss_outcome(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "H",cc00,cc01,cc10,cc11,rho)


pdf(file=glue::glue("plots/BART_LReg_f2_f2com_diff_forms/simulation_loss_outcome_LReg.pdf"),width = 16,height = 15,onefile=FALSE) 
print(ggarrange(
  p.A[[1]], p.B[[1]], p.C[[1]], p.D[[1]],
  p.E[[1]], p.F[[1]], p.G[[1]], p.H[[1]],
  p.A[[2]], p.B[[2]], p.C[[2]], p.D[[2]],
  p.E[[2]], p.F[[2]], p.G[[2]], p.H[[2]],
  labels=c(
    'A','B','C','D',
    'E','F','G','H',
    'A','B','C','D',
    'E','F','G','H'),
  nrow = 4,ncol=4,common.legend = TRUE,legend = 'bottom'))
dev.off() 


#===============================
# plot accuracy when threhold from 0~0.04
#===============================
# cc00 <- 1
# cc01 <- 0
# cc10 <- 1
# cc11 <- 0
# rho <- 0
# threhold.max <- 0.4
# step <- 0.05
# p.A <- plot_BART_LReg_acc(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "A",cc00,cc01,cc10,cc11,rho)
# p.B <- plot_BART_LReg_acc(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "B",cc00,cc01,cc10,cc11,rho)
# p.C <- plot_BART_LReg_acc(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "C",cc00,cc01,cc10,cc11,rho)
# p.D <- plot_BART_LReg_acc(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "D",cc00,cc01,cc10,cc11,rho)
# p.E <- plot_BART_LReg_acc(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "E",cc00,cc01,cc10,cc11,rho)
# p.F <- plot_BART_LReg_acc(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "F",cc00,cc01,cc10,cc11,rho)
# p.G <- plot_BART_LReg_acc(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "G",cc00,cc01,cc10,cc11,rho)
# p.H <- plot_BART_LReg_acc(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "H",cc00,cc01,cc10,cc11,rho)
# 
# pdf(file=glue::glue("plots/BART_LReg_f2_f2com_diff_forms/simulation_acc_LReg.pdf"),width = 15,height = 8,onefile=FALSE) 
# print(ggarrange(
#   p.A, p.B, p.C, p.D,
#   p.E, p.F, p.G, p.H,
#   labels=c(
#     'A','B','C','D',
#     'E','F','G','H'),
#   nrow = 2,ncol=4,common.legend = TRUE,legend = 'bottom'))
# dev.off() 
# 

# =================================================
# plot accuracy, precision, recall for threhold from 0 -0.4
# =================================================
filename <- 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_'
threhold.max <- 0.4
step <- 0.05
ct00 <- 1
ct01 <- 1
ct10 <- 0
ct11 <- 0
cc00 <- 1
cc01 <- 0
cc10 <- 1
cc11 <- 0
rho <- 0

p.A <- plot_BART_LReg_acc_precision_recall(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "A",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
p.B <- plot_BART_LReg_acc_precision_recall(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "B",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
p.C <- plot_BART_LReg_acc_precision_recall(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "C",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
p.D <- plot_BART_LReg_acc_precision_recall(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "D",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
p.E <- plot_BART_LReg_acc_precision_recall(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "E",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
p.F <- plot_BART_LReg_acc_precision_recall(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "F",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
p.G <- plot_BART_LReg_acc_precision_recall(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "G",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
p.H <- plot_BART_LReg_acc_precision_recall(threhold.max = threhold.max, step = 0.05, filename = 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_', scenario = "H",ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)

pdf(file=glue::glue("plots/simulation_BART_LReg_acc_precision_recall.pdf"),width = 15,height = 20,onefile=FALSE) 
print(ggarrange(
  p.A[[1]], p.B[[1]], p.C[[1]], p.D[[1]],
  p.E[[1]], p.F[[1]], p.G[[1]], p.H[[1]],
  p.A[[2]], p.B[[2]], p.C[[2]], p.D[[2]],
  p.E[[2]], p.F[[2]], p.G[[2]], p.H[[2]],
  p.A[[3]], p.B[[3]], p.C[[3]], p.D[[3]],
  p.E[[3]], p.F[[3]], p.G[[3]], p.H[[3]],
  labels=c(
    'A','B','C','D',
    'E','F','G','H',
    'A','B','C','D',
    'E','F','G','H',
    'A','B','C','D',
    'E','F','G','H'),
  nrow = 6,ncol=4,common.legend = TRUE,legend = 'bottom'))
dev.off() 


# =================================================
# print precision,recall,accuracy,loss, value when threhold=0
# =================================================
# results.precision <- results.recall <- results.acc <- loss <- results.outcome <- data.frame(matrix(0,nrow = 0, ncol =5 ))
# colnames(results.precision) <- colnames(results.recall) <- colnames(results.acc) <- colnames(loss) <- colnames(results.outcome) <-c("weight","scenario","model","mean","sd")
# Results.acc <- Results.precision <- Results.recall <- data.frame(matrix(0,nrow = 0, ncol =5))
# 
# filename <- 'results/BART_LReg_f2_f2com_diff_forms/BART_LReg_'
# threhold.max <- 0.4
# step <- 0.05
# ct00 <- 1+threhold
# ct01 <- 1+threhold
# ct10 <- threhold
# ct11 <- threhold
# cc00 <- 1
# cc01 <- 0
# cc10 <- 1
# cc11 <- 0
# rho <- 0
# 
# i <- 0
# 
# for(threhold in seq(0,threhold.max,step)){
#   for (scenario in c("A","B","C","D","E","F","G","H")) {
#     # load the data, which has a list called "outcome"
#     load(glue::glue(filename,threhold,scenario,'.RData'))
#     
#     results.scenario <- print_acc_metrics_BARTLReg(outcome,
#                                                    ct00,ct01,ct10,ct11,
#                                                    cc00,cc01,cc10,cc11,rho)
#     results.scenario.acc <- get_acc(outcome,threhold,cc00,cc01,cc10,cc11,rho)
#     loss_outcome <- print_loss_outcome_BARTLReg(outcome,ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
#     
#     i <- i+1
#     
#     results.precision[i,"weight"] <- threhold
#     results.precision[i,"scenario"] <- scenario
#     results.precision[i,"model"] <- "f1"
#     results.precision[i,"mean"] <- results.scenario[[1]][1]
#     results.precision[i,"sd"] <- results.scenario[[1]][4]
#     results.recall[i,"weight"] <- threhold
#     results.recall[i,"scenario"] <- scenario
#     results.recall[i,"model"] <- "f1"
#     results.recall[i,"mean"] <- results.scenario[[4]][1]
#     results.recall[i,"sd"] <- results.scenario[[4]][4]
#     results.acc[i,"weight"] <- threhold
#     results.acc[i,"scenario"] <- scenario
#     results.acc[i,"model"] <- "f1"
#     results.acc[i,"mean"] <- results.scenario.acc[[1]][2]
#     results.acc[i,"sd"] <- results.scenario.acc[[1]][4]
#     loss[i,"weight"] <- threhold
#     loss[i,"scenario"] <- scenario
#     loss[i,"model"] <- "f1"
#     loss[i,"mean"] <- loss_outcome[[1]][1]
#     loss[i,"sd"] <- loss_outcome[[1]][4]
#     results.outcome[i,"weight"] <- threhold
#     results.outcome[i,"scenario"] <- scenario
#     results.outcome[i,"model"] <- "f1"
#     results.outcome[i,"mean"] <- loss_outcome[[4]][1]
#     results.outcome[i,"sd"] <- loss_outcome[[4]][4]
#     
#     i <- i +1
#     results.precision[i,"weight"] <- threhold
#     results.precision[i,"scenario"] <- scenario
#     results.precision[i,"model"] <- "f2"
#     results.precision[i,"mean"] <- results.scenario[[2]][1]
#     results.precision[i,"sd"] <- results.scenario[[2]][4]
#     results.recall[i,"weight"] <- threhold
#     results.recall[i,"scenario"] <- scenario
#     results.recall[i,"model"] <- "f2"
#     results.recall[i,"mean"] <- results.scenario[[5]][1]
#     results.recall[i,"sd"] <- results.scenario[[5]][4]
#     results.acc[i,"weight"] <- threhold
#     results.acc[i,"scenario"] <- scenario
#     results.acc[i,"model"] <- "f2"
#     results.acc[i,"mean"] <- results.scenario.acc[[2]][2]
#     results.acc[i,"sd"] <- results.scenario.acc[[2]][4]
#     loss[i,"weight"] <- threhold
#     loss[i,"scenario"] <- scenario
#     loss[i,"model"] <- "f2"
#     loss[i,"mean"] <- loss_outcome[[2]][1]
#     loss[i,"sd"] <- loss_outcome[[2]][4]
#     results.outcome[i,"weight"] <- threhold
#     results.outcome[i,"scenario"] <- scenario
#     results.outcome[i,"model"] <- "f2"
#     results.outcome[i,"mean"] <- loss_outcome[[5]][1]
#     results.outcome[i,"sd"] <- loss_outcome[[5]][4]
#     
#     i <- i +1
#     results.precision[i,"weight"] <- threhold
#     results.precision[i,"scenario"] <- scenario
#     results.precision[i,"model"] <- "f3"
#     results.precision[i,"mean"] <- results.scenario[[3]][1]
#     results.precision[i,"sd"] <- results.scenario[[3]][4]
#     results.recall[i,"weight"] <- threhold
#     results.recall[i,"scenario"] <- scenario
#     results.recall[i,"model"] <- "f3"
#     results.recall[i,"mean"] <- results.scenario[[6]][1]
#     results.recall[i,"sd"] <- results.scenario[[6]][4]
#     results.acc[i,"weight"] <- threhold
#     results.acc[i,"scenario"] <- scenario
#     results.acc[i,"model"] <- "f3"
#     results.acc[i,"mean"] <- results.scenario.acc[[3]][2]
#     results.acc[i,"sd"] <- results.scenario.acc[[3]][4]
#     loss[i,"weight"] <- threhold
#     loss[i,"scenario"] <- scenario
#     loss[i,"model"] <- "f3"
#     loss[i,"mean"] <- loss_outcome[[3]][1]
#     loss[i,"sd"] <- loss_outcome[[3]][4]
#     results.outcome[i,"weight"] <- threhold
#     results.outcome[i,"scenario"] <- scenario
#     results.outcome[i,"model"] <- "f3"
#     results.outcome[i,"mean"] <- loss_outcome[[6]][1]
#     results.outcome[i,"sd"] <- loss_outcome[[6]][4]
#   }
#   
#   Results.acc <- rbind(Results.acc, results.acc)
#   Results.precision <- rbind(Results.precision, results.precision)
#   Results.recall <- rbind(Results.recall, results.recall)
# }
# 

# results.acc <- cbind(results.acc,acc.mean.upper.lower)[,c("scenario","model","acc.mean.upper.lower")]
# 
# precision.mean.upper.lower <- apply(results.precision, 1, FUN=function(x) paste(round(as.numeric(x[3]),digits = 3),"(",round(as.numeric(x[4]),digits = 3),",",round(as.numeric(x[5]),digits = 3),")"))
# results.precision <- cbind(results.precision,precision.mean.upper.lower)[,c("scenario","model","precision.mean.upper.lower")]
# 
# recall.mean.upper.lower <- apply(results.recall, 1, FUN=function(x) paste(round(as.numeric(x[3]),digits = 3),"(",round(as.numeric(x[4]),digits = 3),",",round(as.numeric(x[5]),digits = 3),")"))
# results.recall <- cbind(results.recall,recall.mean.upper.lower)[,c("scenario","model","recall.mean.upper.lower")]
# 
# loss.mean.upper.lower <- apply(loss, 1, FUN=function(x) paste(round(as.numeric(x[3]),digits = 3),"(",round(as.numeric(x[4]),digits = 3),",",round(as.numeric(x[5]),digits = 3),")"))
# loss <- cbind(loss,loss.mean.upper.lower)[,c("scenario","model","loss.mean.upper.lower")]
# 
# outcome.mean.upper.lower <- apply(results.outcome, 1, FUN=function(x) paste(round(as.numeric(x[3]),digits = 3),"(",round(as.numeric(x[4]),digits = 3),",",round(as.numeric(x[5]),digits = 3),")"))
# results.outcome <- cbind(results.outcome,outcome.mean.upper.lower)[,c("scenario","model","outcome.mean.upper.lower")]
# 
# cbind(loss$scenario,"&",loss$model,"&",results.outcome$outcome.mean.upper.lower,"&",results.acc$acc.mean.upper.lower,"&",results.precision$precision.mean.upper.lower,"&",results.recall$recall.mean.upper.lower, "\\")
# 
# cbind(loss,results.outcome$outcome.mean.upper.lower,results.acc$acc.mean.upper.lower,results.precision$precision.mean.upper.lower,results.recall$recall.mean.upper.lower)
# 
# threhold <- 0
# depth <- 1
# filename <- 'results/BART_TREE10000_'
# load(paste(filename,'A.RData',sep='')); loss.outcome.A <- get_BART_TREE_loss_outcome(threhold, depth, outcome,scenario="A")
# load(paste(filename,'B.RData',sep='')); loss.outcome.B <- get_BART_TREE_loss_outcome(threhold, depth, outcome,scenario="B")
# load(paste(filename,'C.RData',sep='')); loss.outcome.C <- get_BART_TREE_loss_outcome(threhold, depth, outcome,scenario="C")
# load(paste(filename,'D.RData',sep='')); loss.outcome.D <- get_BART_TREE_loss_outcome(threhold, depth, outcome,scenario="D")
# load(paste(filename,'E.RData',sep='')); loss.outcome.E <- get_BART_TREE_loss_outcome(threhold, depth, outcome,scenario="E")
# load(paste(filename,'F.RData',sep='')); loss.outcome.F <- get_BART_TREE_loss_outcome(threhold, depth, outcome,scenario="F")
# load(paste(filename,'G.RData',sep='')); loss.outcome.G <- get_BART_TREE_loss_outcome(threhold, depth, outcome,scenario="G")
# load(paste(filename,'H.RData',sep='')); loss.outcome.H <- get_BART_TREE_loss_outcome(threhold, depth, outcome,scenario="H")
# loss.tree <- rbind(loss.outcome.A$loss,
#                    loss.outcome.B$loss,
#                    loss.outcome.C$loss,
#                    loss.outcome.D$loss,
#                    loss.outcome.E$loss,
#                    loss.outcome.F$loss,
#                    loss.outcome.G$loss,
#                    loss.outcome.H$loss)





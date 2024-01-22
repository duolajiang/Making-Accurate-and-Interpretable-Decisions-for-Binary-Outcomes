DataImportGeneral <- function(year){
  data <- data_process_all_variable_update()
  data <- data_process_to_dummy_update(data)
  data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
  data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
  # target population index
  if(year==5){
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&
      ((data$Stage==2)|(data$Stage==3))&
      (data$year_diagnosis<2014)&
      (data$vitfup>0.082)
    data <- data[target_index,]
  } else if (year==2){
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&
      ((data$Stage==2)|(data$Stage==3)) &
      (data$vitfup>0.082)
    data <- data[target_index,]
  } else{ #for 20156,follow up year=2
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&
      ((data$Stage==2)|(data$Stage==3))&
      (data$year_diagnosis>2014) &
      (data$vitfup>0.082)
    data <- data[target_index,]
  }
  return(data)
}

# =============================================================
# compute posterior expected loss and outcome based on posterior distribution of theta
# =============================================================
posterior_computation <- function(data,xvars,year,ntree){
  ## there is not much difference in estiamting Y=1(death) or Y=1(alive)
  set.seed(123)
  if((year==2)|(year==5)){
    length_test <- dim(data)[1]
    train <- data[,c(xvars)]
    y <- data$combined_chemo
    test <- data[,xvars]
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    data$bartps <- apply(post$prob.test,2,mean)

    selection <- data
    selection <- filter(selection,((vitfup<year)&(vitstat==1))|(vitfup>=year))
    selection[selection$vitfup>=year,"vitstat"] <- 0

    length_test <- dim(data)[1]
    train <- selection[,c('combined_chemo',xvars,'bartps')]
    y <- selection$vitstat
    test <- data[,c('combined_chemo',xvars,'bartps')]
    test_bart_1 <- test_bart_0 <- test
    test_bart_1$combined_chemo <- 1
    test_bart_0$combined_chemo <- 0
    test <- rbind(test_bart_1,test_bart_0)
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    post$prob.test <- 1- post$prob.test
    post_y1_y0 <- post$prob.test
    file.name <- ifelse(year==5,'post_y1_y0_5fp.RData','post_y1_y0_2fp.RData')
    save(post_y1_y0,file=file.name)
  } else{
    length_test <- dim(data)[1]
    train <- data[,c(xvars)]
    y <- data$combined_chemo
    test <- data[,xvars]
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    ps <- apply(post$prob.test,2,mean)
    data$bartps <- ps

    year <- 2
    selection <- data
    selection <- filter(selection,((vitfup<year)&(vitstat==1))|(vitfup>=year))
    selection[selection$vitfup>=year,"vitstat"] <- 0

    length_test <- dim(data)[1]
    train <- selection[,c('combined_chemo',xvars,'bartps')]
    y <- selection$vitstat
    test <- data[,c('combined_chemo',xvars,'bartps')]
    test_bart_1 <- test_bart_0 <- test
    test_bart_1$combined_chemo <- 1
    test_bart_0$combined_chemo <- 0
    test <- rbind(test_bart_1,test_bart_0)
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    post$prob.test <- 1- post$prob.test
    post_y1_y0 <- post$prob.test
    save(post_y1_y0, file = "post_y1_y0_2fp_20156.RData")
  }
}





# =============================================================
# xvars.all: the covariates used to estimate the potential outcomes setup when load_data=FALSE
# year: 5fp or 2fp analysis, if year=5, then only data before 2014 then is used
# load_data: TRUE load existing potential outcomes or FALSE computing new potential outcomes
# alpha: significance level for splitting criteria
# min.ndsz: minimum number of observations in each node
# n0: min(nl.t,nl.c,nr.t,nr.c)
# max.depth: max depth of tree
# rho: the correlation of potential outcomes
# cd: cost for death
# ct_seq: cost for treatment, either ct_seq= seq(0,1,0.01) or fixed value, like 0.05
# =============================================================
PolicyEvaluation <- function(xvars.all=NULL, year, load_data, alpha=0.05, min.ndsz=200, n0=100, max.depth=15, rho=0, cd=1, ct_seq=0.05){
  # =============================================================
  # PREPARATION: IMPORT DATA, COMPUTE GUIDELINE, SELECT TARGET POPULATION, AND POTENTIAL OUTCOME COMPUTATIONS
  # =============================================================
  ntree <- 50
  data <- data_process_all_variable()
  data <- data_process_to_dummy(data)
  data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
  data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
  # target population index
  if(year==5){
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis<2014)
    data <- data[target_index,]
  } else{
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))
    data <- data[target_index,]
  }

  # load/compute potential outcomes
  if(load_data==TRUE){
    file.name <- ifelse(year==5,'post_y1_y0_5fp.RData','post_y1_y0_2fp.RData')
    load(file.name)
  } else{
    post_y1_y0 <- posterior_computation(data,xvars.all,year,ntree) #nrow=1000, ncol=dim(data)*2
    file.name <- ifelse(year==5,'post_y1_y0_5fp.RData','post_y1_y0_2fp.RData')
    save(post_y1_y0,file=file.name)
  }


  length_test <- (dim(post_y1_y0)[2])/2
  ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
  ate.test.mean <- apply(ate.test,2,mean)
  p1_ <- post_y1_y0[,c(1:length_test)]
  p_1 <- post_y1_y0[,c((length_test+1):(2*length_test))]

  p1_ <- p1_[,c(target_index)]
  p_1 <- p_1[,c(target_index)]
  ate.test.mean <- ate.test.mean[c(target_index)]
  expected.y1 <- apply(p1_, 2, mean)
  expected.y0 <- apply(p_1, 2, mean)

  #target population adding ate for each observation
  data <- data[c(target_index),]
  #combine lymphomatic invation and agio invation together
  data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                           ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
  data$ate <- ate.test.mean
  data$y1 <- expected.y1
  data$y0 <- expected.y0

  # =============================================================
  # POLICY FOR INTERPRETABLE TREE
  # =============================================================
  xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','num_cormor','pT_num',
              'Stage','sex','MSI','pN_num','cM_num',
              "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc'
  )
  dat <- data[,xvars]
  dat$trt <- data$combined_chemo
  dat$y <- ate.test.mean
  # COLUMNS OF COVARIATES
  split.var <- 1:16;ctg=7:16; mtry=length(split.var)
  test=NULL;
  # FLEXIBILTIY OF TREE STRUCTURE
  min.ndsz=min.ndsz # mininum number of observations in a node
  n0=n0 # minimun number of treatment/control in a node
  max.depth=max.depth # max depth of tree
  alpha=alpha # only when significance level larger than 0.05, then the split is taken

  tree <- grow.INT(dat, test=NULL,
                   split.var,ctg,
                   min.ndsz, n0, max.depth,alpha,
                   mtry=length(split.var))

  tree.ate.hat <- predict_initree_ITE(tree,dat,ctg)

  # =============================================================
  # correlation of Y(1),Y(0): rho
  # theta
  # =============================================================
  rho <- rho
  p11 <- rho*sqrt(p1_*(1-p1_)*p_1*(1-p_1)) + p1_*p_1
  p10 <- p1_ - p11
  p01 <- p_1 - p11
  p00 <- 1- p11 - p10 - p01
  # =============================================================
  # COMPUTE POSTERIOR EXPECTED LOSS AND POSTERIOR EXPECTED OUTCOME
  # =============================================================
  cd <- cd
  ct_seq <- ct_seq
  guideline2014.action <- ifelse(data$guideline2014=='unknown',ifelse(data$ate>0,1,0),data$guideline2014)
  guideline2018.action <- ifelse(data$guideline2018=='unknown',ifelse(data$ate>0,1,0),data$guideline2018)
  observation.action <- data$combined_chemo


  filename <- ifelse(year==5,ifelse(length(ct_seq)==1, paste('loss_outcome_value_5fp_',ct_seq,'.RData',sep=''),'loss_outcome_value_all_5fp.RData'),
                     ifelse(length(ct_seq)==1, paste('loss_outcome_value_2fp_',ct_seq,'.RData',sep=''),'loss_outcome_value_all_2fp.RData'))

  expected_loss_outcome_estimation(rho,p11,p10,p01,p00,cd,ct_seq,data$ate,tree.ate.hat,
                                   guideline2014.action,guideline2018.action,observation.action,expected.y1,expected.y0,
                                   filename=filename)

}

PolicyEvaluationUpdate <- function(xvars.all=NULL, year, load_data, alpha=0.05, min.ndsz=200, n0=100, max.depth=15, rho=0, cd=1, ct_seq=0.05){
  # =============================================================
  # PREPARATION: IMPORT DATA, COMPUTE GUIDELINE, SELECT TARGET POPULATION, AND POTENTIAL OUTCOME COMPUTATIONS
  # =============================================================
  #browser()
  ntree <- 50
  data <- data_process_all_variable_update()
  data <- data_process_to_dummy_update(data)
  data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
  data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
  # target population index
  if(year==5){
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis<2014)
    data <- data[target_index,]
  } else{
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))
    data <- data[target_index,]
  }

  # load/compute potential outcomes
  if(load_data==TRUE){
    file.name <- ifelse(year==5,'post_y1_y0_5fp_update.RData','post_y1_y0_2fp_update.RData')
    load(file.name)
  } else{
    post_y1_y0 <- posterior_computation(data,xvars.all,year,ntree) #nrow=1000, ncol=dim(data)*2
    file.name <- ifelse(year==5,'post_y1_y0_5fp_update.RData','post_y1_y0_2fp_update.RData')
    save(post_y1_y0,file=file.name)
  }


  length_test <- (dim(post_y1_y0)[2])/2
  ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
  ate.test.mean <- apply(ate.test,2,mean)
  p1_ <- post_y1_y0[,c(1:length_test)]
  p_1 <- post_y1_y0[,c((length_test+1):(2*length_test))]
  expected.y1 <- apply(p1_, 2, mean)
  expected.y0 <- apply(p_1, 2, mean)

  #combine lymphomatic invation and agio invation together
  data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                           ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
  data$ate <- ate.test.mean
  data$y1 <- expected.y1
  data$y0 <- expected.y0

  #=============================================================
  #POLICY FOR INTERPRETABLE TREE
  #=============================================================
  xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num','year_diagnosis',
              'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
              "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
  )
  dat <- data[,xvars]
  dat$trt <- data$combined_chemo
  dat$y <- ate.test.mean
  # COLUMNS OF COVARIATES
  split.var <- 1:20;ctg=7:20; mtry=length(split.var)
  test=NULL;
  # FLEXIBILTIY OF TREE STRUCTURE
  min.ndsz=min.ndsz # mininum number of observations in a node
  n0=n0 # minimun number of treatment/control in a node
  max.depth=max.depth # max depth of tree
  alpha=alpha # only when significance level larger than 0.05, then the split is taken

  tree <- grow.INT(dat, test=NULL,
                   split.var,ctg,
                   min.ndsz, n0, max.depth,alpha,
                   mtry=length(split.var))

  tree.ate.hat <- predict_initree_ITE(tree,dat,ctg)

  # =============================================================
  # correlation of Y(1),Y(0): rho
  # theta
  # =============================================================
  rho <- rho
  p11 <- rho*sqrt(p1_*(1-p1_)*p_1*(1-p_1)) + p1_*p_1
  p10 <- p1_ - p11
  p01 <- p_1 - p11
  p00 <- 1- p11 - p10 - p01
  # =============================================================
  # COMPUTE POSTERIOR EXPECTED LOSS AND POSTERIOR EXPECTED OUTCOME
  # =============================================================
  cd <- cd
  ct_seq <- ct_seq
  guideline2014.action <- ifelse(data$guideline2014=='unknown',ifelse(data$ate>0,1,0),data$guideline2014)
  guideline2018.action <- ifelse(data$guideline2018=='unknown',ifelse(data$ate>0,1,0),data$guideline2018)
  observation.action <- data$combined_chemo


  filename <- ifelse(year==5,ifelse(length(ct_seq)==1, paste('loss_outcome_value_5fp_update',ct_seq,'.RData',sep=''),'loss_outcome_value_all_update_5fp.RData'),
                     ifelse(length(ct_seq)==1, paste('loss_outcome_value_2fp_update',ct_seq,'.RData',sep=''),'loss_outcome_value_all_update_2fp.RData'))

  expected_loss_outcome_estimation(rho,p11,p10,p01,p00,cd,ct_seq,data$ate,tree.ate.hat,
                                   guideline2014.action,guideline2018.action,observation.action,expected.y1,expected.y0,
                                   filename=filename)

}

TreatmentProprotionAndConfusionMatrix <- function(cost,year,alpha=0.05, min.ndsz=200, n0=100, max.depth=15){
  ntree <- 50
  data <- data_process_all_variable_update()
  data <- data_process_to_dummy_update(data)
  data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
  data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
  # target population index
  if(year==5){
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis<2014)
    data <- data[target_index,]
  } else{
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))
    data <- data[target_index,]
  }

  file.name <- ifelse(year==5,'post_y1_y0_5fp_update.RData','post_y1_y0_2fp_update.RData')
  load(file.name)

  length_test <- (dim(post_y1_y0)[2])/2
  ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
  ate.test.mean <- apply(ate.test,2,mean)

  #combine lymphomatic invation and agio invation together
  data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                           ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
  data$ate <- ate.test.mean
  #=============================================================
  #POLICY FOR INTERPRETABLE TREE
  #=============================================================
  xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num','year_diagnosis',
              'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
              "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
  )
  dat <- data[,xvars]
  dat$trt <- data$combined_chemo
  dat$y <- ate.test.mean
  # COLUMNS OF COVARIATES
  split.var <- 1:20;ctg=7:20; mtry=length(split.var)
  test=NULL;
  # FLEXIBILTIY OF TREE STRUCTURE
  min.ndsz=min.ndsz # mininum number of observations in a node
  n0=n0 # minimun number of treatment/control in a node
  max.depth=max.depth # max depth of tree
  alpha=alpha # only when significance level larger than 0.05, then the split is taken

  tree <- grow.INT(dat, test=NULL,
                   split.var,ctg,
                   min.ndsz, n0, max.depth,alpha,
                   mtry=length(split.var))

  tree.ate.hat <- predict_initree_ITE(tree,dat,ctg)

  tree.action <- round(mean(ifelse(tree.ate.hat>cost,1,0))*100,digits = 1)
  optimal.action <- round(mean(ifelse(dat$y>cost,1,0))*100,digits = 1)
  guideline2014.action <- round(mean(as.numeric(ifelse(data$guideline2014=='unknown',ifelse(data$ate>0,1,0),data$guideline2014))),digits = 1)
  guideline2018.action <- round(mean(as.numeric(ifelse(data$guideline2018=='unknown',ifelse(data$ate>0,1,0),data$guideline2018))),digits = 1)
  observation.action <- round(mean(data$combined_chemo)*100,digits = 1)
  return(cbind(c('tree','optimal','guideline2018','guideline2014','observation'),
               c(tree.action,optimal.action,guideline2018.action,guideline2014.action,observation.action)))
}
# =============================================================
# compute and visualize the posterior expected loss, outcome, value for each individual with 95%credible interval
# as cost=0.05,0.1,...
# =============================================================
PolicyUncertaintyEvaluation <- function(xvars.all=NULL, year, load_data, alpha=0.05, min.ndsz=200, n0=100, max.depth=15, rho=0, cd=1, ct=0.05){
  # =============================================================
  # PREPARATION: IMPORT DATA, COMPUTE GUIDELINE, SELECT TARGET POPULATION, AND POTENTIAL OUTCOME COMPUTATIONS
  # =============================================================
  #browser()
  ntree <- 50
  data <- data_process_all_variable_update()
  data <- data_process_to_dummy_update(data)
  data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
  data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
  # target population index
  if(year==5){
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis<2014)
    data <- data[target_index,]
  } else{
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))
    data <- data[target_index,]
  }

  # load/compute potential outcomes
  file.name <- ifelse(year==5,'post_y1_y0_5fp_update.RData','post_y1_y0_2fp_update.RData')
  load(file.name)

  length_test <- (dim(post_y1_y0)[2])/2
  ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
  ate.test.mean <- apply(ate.test,2,mean)
  p1_ <- post_y1_y0[,c(1:length_test)]
  p_1 <- post_y1_y0[,c((length_test+1):(2*length_test))]
  expected.y1 <- apply(p1_, 2, mean)
  expected.y0 <- apply(p_1, 2, mean)
  data$ate <- ate.test.mean

  #combine lymphomatic invation and agio invation together
  data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                           ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))

  #=============================================================
  #POLICY FOR INTERPRETABLE TREE
  #=============================================================
  xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num','year_diagnosis',
              'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
              "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
  )
  dat <- data[,xvars]
  dat$trt <- data$combined_chemo
  dat$y <- ate.test.mean
  # COLUMNS OF COVARIATES
  split.var <- 1:20;ctg=7:20; mtry=length(split.var)
  test=NULL;
  # FLEXIBILTIY OF TREE STRUCTURE
  min.ndsz=min.ndsz # mininum number of observations in a node
  n0=n0 # minimun number of treatment/control in a node
  max.depth=max.depth # max depth of tree
  alpha=alpha # only when significance level larger than 0.05, then the split is taken

  tree <- grow.INT(dat, test=NULL,
                   split.var,ctg,
                   min.ndsz, n0, max.depth,alpha,
                   mtry=length(split.var))

  tree.ate.hat <- predict_initree_ITE(tree,dat,ctg)

  # =============================================================
  # correlation of Y(1),Y(0): rho
  # theta00, theta01, theta10, theta11
  # =============================================================
  rho <- rho
  p11 <- rho*sqrt(p1_*(1-p1_)*p_1*(1-p_1)) + p1_*p_1
  p10 <- p1_ - p11
  p01 <- p_1 - p11
  p00 <- 1- p11 - p10 - p01
  # =============================================================
  # COMPUTE POSTERIOR EXPECTED LOSS AND POSTERIOR EXPECTED OUTCOME
  # =============================================================
  guideline2014.action <- as.numeric(ifelse(data$guideline2014=='unknown',ifelse(data$ate>0,1,0),data$guideline2014))
  guideline2018.action <- as.numeric(ifelse(data$guideline2018=='unknown',ifelse(data$ate>0,1,0),data$guideline2018))
  observation.action <- data$combined_chemo

  output <- expected_loss_outcome_value_uncertainty(rho,p1_,p_1,cd=cd,ct=ct,tree.ate.hat,guideline2018.action,observation.action)
  return(output)
}

# =============================================================
# cost-sensitive tree visualization for a) old data b) updated data c) 2015&2016
# =============================================================
# year: 5fp or 2fp analysis, if year=5, then only data before 2014 then is used
# alpha: significance level for splitting criteria
# min.ndsz: minimum number of observations in each node
# n0: min(nl.t,nl.c,nr.t,nr.c)
# max.depth: max depth of tree
# return tree, dat, and cate categorical covariates
# the single tree computation, the same as tree computation in PolicyEvaluation
# =============================================================
DataImportAndTreeGrown <- function(year,alpha=0.05,min.ndsz=200, pt=pt, max.depth=15){
  data <- data_process_all_variable()
  data <- data_process_to_dummy(data)
  # target population index
  if(year==5){
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis<2014)
    load('post_y1_y0_5fp.RData')
  } else{
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))
    load('post_y1_y0_2fp.RData')
  }

  length_test <- (dim(post_y1_y0)[2])/2
  ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
  ate.test.mean <- apply(ate.test,2,mean)
  ate.test.mean <- ate.test.mean[c(target_index)]
  data <- data[c(target_index),]
  data$lym_agioinvasie <- ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                 ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other'))
  data$ate <- ate.test.mean

  # xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','num_cormor','pT_num',
  #             'Stage','sex','MSI','pN_num','cM_num',
  #             "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc'
  # )
  xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num',
              'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
              "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
  )
  dat <- data[,xvars]
  dat$trt <- data$combined_chemo
  dat$y <- ate.test.mean
  # COLUMNS OF COVARIATES
  split.var <- 1:19;ctg=6:19; mtry=length(split.var)
  #split.var <- 1:16;ctg=7:16; mtry=length(split.var)
  test=NULL;
  # FLEXIBILTIY OF TREE STRUCTURE
  min.ndsz=min.ndsz # mininum number of observations in a node
  pt=pt # minimun number of treatment/control in a node
  max.depth=max.depth # max depth of tree
  alpha=alpha # only when significance level larger than 0.05, then the split is taken

  tree <- grow.INT(dat, test=NULL,
                   split.var,ctg,
                   min.ndsz, pt=pt, max.depth,alpha,
                   mtry=length(split.var))

  varname <- as.character(unique(tree$vname[!is.na(tree$vname)]))
  cate <- NULL
  for (var in varname){
    if(class(dat[,var])=='factor'){
      cate <- cbind(cate,var)
    }
  }

  return(list(tree,cate,dat))
}

DataImportAndTreeGrownUpdate <- function(year,stage=NULL,alpha=0.05,min.ndsz=200, pt=pt, max.depth=15){
  data <- data_process_all_variable_update()
  data <- data_process_to_dummy_update(data)

  # target population index
  if(year==5){
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis<2014)
    data <- data[c(target_index),]
    load('post_y1_y0_5fp_update.RData')
  } else{
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))
    data <- data[c(target_index),]
    load('post_y1_y0_2fp_update.RData')
  }


  length_test <- (dim(post_y1_y0)[2])/2
  ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
  ate.test.mean <- apply(ate.test,2,mean)
  data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                           ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
  data$ate <- ate.test.mean
  data$predictive_var <- apply(ate.test,2,var)

  if(!is.null(stage)){
    data <- filter(data,Stage==stage)
    xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num',
                'sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
                "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
    )
  } else{
    xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num',
                'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
                "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
    )
  }

  dat <- data[,xvars]
  dat$trt <- data$combined_chemo
  dat$y <- data$ate
  # COLUMNS OF COVARIATES
  split.var <- 1:length(xvars);ctg=6:length(xvars); mtry=length(split.var)
  test=NULL;
  # FLEXIBILTIY OF TREE STRUCTURE
  min.ndsz=min.ndsz # mininum number of observations in a node
  pt=pt # minimun number of treatment/control in a node
  max.depth=max.depth # max depth of tree
  alpha=alpha # only when significance level larger than 0.05, then the split is taken

  tree <- grow.INT(dat, test=NULL,
                   split.var,ctg,
                   min.ndsz, pt=pt, max.depth,alpha,
                   mtry=length(split.var))

  varname <- as.character(unique(tree$vname[!is.na(tree$vname)]))
  cate <- NULL
  for (var in varname){
    if(class(dat[,var])=='factor'){
      cate <- cbind(cate,var)
    }
  }

  return(list(tree,cate,dat))
}
# =============================================================
# evidence of 2fp for 2015 and 2016 data only
# =============================================================
Data20152016ImportAndTreeGrown <- function(alpha=0.05,min.ndsz=200, pt=pt, max.depth=15){
  # =============================================================
  # PREPARATION: IMPORT DATA, COMPUTE GUIDELINE, SELECT TARGET POPULATION, AND POTENTIAL OUTCOME COMPUTATIONS
  # =============================================================
  set.seed(123)
  xvars.all  <- c('age_at_diagnosis','lymph_assessed',
                  'male',#'female',
                  'BMI',
                  'MSI_0','MSI_1','MSI_unknown',#'MSI_missing',
                  'colon_perforation_0',"colon_perforation_1","colon_perforation_unknown",
                  #"colon_perforation_missing",
                  'lymphomatic_invation_no',
                  'lymphomatic_invation_yes',
                  'lymphomatic_invation_suspect',
                  'lymphomatic_invation_unknown',
                  #'lymphatic_vessel_invasion_missing'
                  'agio_invation_no',
                  'agio_invation_EMVI',
                  'agio_invation_IMVI',
                  'agio_invation_suspect',
                  #'agio_invation_NA',
                  'agio_invation_unknown',
                  "Grade_good_moderate","Grade_poor_undifferentiated",
                  #'pT0',
                  'pT1','pT2','pT3',#'pT4',
                  'pN0','pN1','pN2','pNX',#'pN.missing',
                  'cM0','cM1','cMX',#'cM.missing',
                  'max_LoS',
                  'Morphology_1','Morphology_2','Morphology_3','Morphology_4',
                  'Stage2','Stage3',
                  'highrisk_false','highrisk_true','highrisk_unknown', #Morphology_5
                  'proximal','distal',
                  'BRAF_0','BRAF_1','BRAF_9',
                  'RAS_0','RAS_1','RAS_9',
                  'ASAclass_1','ASAclass_2','ASAclass_3','ASAclass_4','ASAclass_5','ASAclass_9',
                  'PS_notregistered','PS_unknown','PS_WHO1','PS_WHO2','PS_WHO3','PS_WHO4'
  )
  if(file_test("-f", "post_y1_p0_2fp_20156.RData")){
    data <- data_process_all_variable_update()
    data <- data_process_to_dummy_update(data)
    data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
    data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis>2014)
    data <- data[target_index,]

    load("post_y1_p0_2fp_20156.RData")
    length_test <- (dim(post_y1_y0)[2])/2
    ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
    ate.test.mean <- apply(ate.test,2,mean)
    #target population adding ate for each observation
    data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                             ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
    data$ate <- ate.test.mean
  } else{
    ntree <- 50
    data <- data_process_all_variable_update()
    data <- data_process_to_dummy_update(data)
    data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
    data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis>2014)
    data <- data[target_index,]

    set.seed(123)
    length_test <- dim(data)[1]
    train <- data[,c(xvars.all)]
    y <- data$combined_chemo
    test <- data[,xvars.all]
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    ps <- apply(post$prob.test,2,mean)
    data$bartps <- ps

    year <- 2
    selection <- data
    selection <- filter(selection,((vitfup<year)&(vitstat==1))|(vitfup>=year))
    selection[selection$vitfup>=year,"vitstat"] <- 0

    length_test <- dim(data)[1]
    train <- selection[,c('combined_chemo',xvars.all,'bartps')]
    y <- selection$vitstat
    test <- data[,c('combined_chemo',xvars.all,'bartps')]
    test_bart_1 <- test_bart_0 <- test
    test_bart_1$combined_chemo <- 1
    test_bart_0$combined_chemo <- 0
    test <- rbind(test_bart_1,test_bart_0)
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    post$prob.test <- 1- post$prob.test

    post_y1_y0 <- post$prob.test
    save(post_y1_y0, file = "post_y1_p0_2fp_20156.RData")
    length_test <- (dim(post_y1_y0)[2])/2
    ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
    ate.test.mean <- apply(ate.test,2,mean)


    #target population adding ate for each observation
    data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                             ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
    data$ate <- ate.test.mean
  }
  # =============================================================
  # POLICY FOR INTERPRETABLE TREE
  # =============================================================
  xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num',
              'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
              "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
  )
  dat <- data[,xvars]
  dat$trt <- data$combined_chemo
  dat$y <- ate.test.mean
  # COLUMNS OF COVARIATES
  split.var <- 1:19;ctg=6:19; mtry=length(split.var)
  test=NULL;
  # FLEXIBILTIY OF TREE STRUCTURE
  min.ndsz=min.ndsz # mininum number of observations in a node
  pt=pt # minimun number of treatment/control in a node
  max.depth=max.depth # max depth of tree
  alpha=alpha # only when significance level larger than 0.05, then the split is taken
  tree <- grow.INT(dat, test=NULL,
                   split.var,ctg,
                   min.ndsz, pt=pt, max.depth,alpha,
                   mtry=length(split.var))


  varname <- as.character(unique(tree$vname[!is.na(tree$vname)]))
  cate <- NULL
  for (var in varname){
    if(class(dat[,var])=='factor'){
      cate <- cbind(cate,var)
    }
  }

  return(list(tree,cate,dat))
}

Data20152016ImportAndTreeGrownILL <- function(min.ndsz=200, pt=0.2, max.depth=15){
  # =============================================================
  # PREPARATION: IMPORT DATA, COMPUTE GUIDELINE, SELECT TARGET POPULATION, AND POTENTIAL OUTCOME COMPUTATIONS
  # =============================================================
  set.seed(123)
  xvars.all  <- c('age_at_diagnosis','lymph_assessed',
                  'male',#'female',
                  'BMI',
                  'MSI_0','MSI_1','MSI_unknown',#'MSI_missing',
                  'colon_perforation_0',"colon_perforation_1","colon_perforation_unknown",
                  #"colon_perforation_missing",
                  'lymphomatic_invation_no',
                  'lymphomatic_invation_yes',
                  'lymphomatic_invation_suspect',
                  'lymphomatic_invation_unknown',
                  #'lymphatic_vessel_invasion_missing'
                  'agio_invation_no',
                  'agio_invation_EMVI',
                  'agio_invation_IMVI',
                  'agio_invation_suspect',
                  #'agio_invation_NA',
                  'agio_invation_unknown',
                  "Grade_good_moderate","Grade_poor_undifferentiated",
                  #'pT0',
                  'pT1','pT2','pT3',#'pT4',
                  'pN0','pN1','pN2','pNX',#'pN.missing',
                  'cM0','cM1','cMX',#'cM.missing',
                  'max_LoS',
                  'Morphology_1','Morphology_2','Morphology_3','Morphology_4',
                  'Stage2','Stage3',
                  'highrisk_false','highrisk_true','highrisk_unknown', #Morphology_5
                  'proximal','distal',
                  'BRAF_0','BRAF_1','BRAF_9',
                  'RAS_0','RAS_1','RAS_9',
                  'ASAclass_1','ASAclass_2','ASAclass_3','ASAclass_4','ASAclass_5','ASAclass_9',
                  'PS_notregistered','PS_unknown','PS_WHO1','PS_WHO2','PS_WHO3','PS_WHO4'
  )
  if(file_test("-f", "post_y1_p0_2fp_20156.RData")){
    data <- data_process_all_variable_update()
    data <- data_process_to_dummy_update(data)
    data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
    data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis>2014)
    data <- data[target_index,]

    load("post_y1_p0_2fp_20156.RData")
    length_test <- (dim(post_y1_y0)[2])/2
    ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
    ate.test.mean <- apply(ate.test,2,mean)
    #target population adding ate for each observation
    data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                             ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
    data$ate <- ate.test.mean
    data$predictive_var <- apply(ate.test,2,var) #var is unbiased estimator of variance
  } else{
    ntree <- 50
    data <- data_process_all_variable_update()
    data <- data_process_to_dummy_update(data)
    data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
    data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis>2014)
    data <- data[target_index,]

    set.seed(123)
    length_test <- dim(data)[1]
    train <- data[,c(xvars.all)]
    y <- data$combined_chemo
    test <- data[,xvars.all]
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    ps <- apply(post$prob.test,2,mean)
    data$bartps <- ps

    year <- 2
    selection <- data
    selection <- filter(selection,((vitfup<year)&(vitstat==1))|(vitfup>=year))
    selection[selection$vitfup>=year,"vitstat"] <- 0

    length_test <- dim(data)[1]
    train <- selection[,c('combined_chemo',xvars.all,'bartps')]
    y <- selection$vitstat
    test <- data[,c('combined_chemo',xvars.all,'bartps')]
    test_bart_1 <- test_bart_0 <- test
    test_bart_1$combined_chemo <- 1
    test_bart_0$combined_chemo <- 0
    test <- rbind(test_bart_1,test_bart_0)
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    post$prob.test <- 1- post$prob.test

    post_y1_y0 <- post$prob.test
    save(post_y1_y0, file = "post_y1_p0_2fp_20156.RData")
    length_test <- (dim(post_y1_y0)[2])/2
    ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
    ate.test.mean <- apply(ate.test,2,mean)


    #target population adding ate for each observation
    data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                             ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
    data$ate <- ate.test.mean
    data$predictive_var <- apply(ate.test,2,var)
  }
  # =============================================================
  # POLICY FOR INTERPRETABLE TREE
  # =============================================================
  xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num',
              'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
              "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
  )
  dat <- data[,xvars]
  dat$trt <- data$combined_chemo
  dat$y <- ate.test.mean
  dat$predictive_var <- data$predictive_var
  # COLUMNS OF COVARIATES
  split.var <- 1:19;ctg=6:19; mtry=length(split.var)
  test=NULL;
  # FLEXIBILTIY OF TREE STRUCTURE
  min.ndsz=min.ndsz # mininum number of observations in a node
  pt=pt # minimun number of treatment/control in a node
  max.depth=max.depth # max depth of tree
  tree <- grow.INT.ILL(dat, test=NULL,
                       split.var,ctg,
                       min.ndsz, pt, max.depth,
                       mtry=length(split.var))


  varname <- as.character(unique(tree$vname[!is.na(tree$vname)]))
  cate <- NULL
  for (var in varname){
    if(class(dat[,var])=='factor'){
      cate <- cbind(cate,var)
    }
  }

  return(list(tree,cate,dat))
}
# =============================================================
# for printing
# =============================================================
GetNodeName <- function(tree,dat,cate,nodeid){
  #browser()
  isleft <- (substr(nodeid,start = nchar(nodeid),stop = nchar(nodeid))==1)
  node.parent.id <- substr(nodeid,start = 1,stop = nchar(nodeid)-1)
  node.parent <- tree[tree$node==node.parent.id,]
  var.name <- as.character(node.parent[,'vname'])
  if(is.element(var.name,cate)){
    cut <- as.character(unlist(strsplit(as.character(node.parent[,'cut']),split=" ")))
    cut.index <- is.element(levels(dat[,var.name]),cut)
    cut.val.left <- levels(dat[,var.name])[cut.index];   cut.val.left <- CombineLevelsString(cut.val.left)
    cut.val.right <- levels(dat[,var.name])[!cut.index]; cut.val.right <- CombineLevelsString(cut.val.right)
    target <- ifelse(isleft==TRUE,paste(var.name,'=',cut.val.left,sep = ''),paste(var.name,'=',cut.val.right,sep = ''))
  } else{
    cut.val <- as.character(node.parent[,'cut'])
    target <- ifelse(isleft==TRUE,paste(var.name,'<=',cut.val,sep = ''),paste(var.name,'>',cut.val,sep = ''))
  }
  return(target)
}
CombineLevelsString <- function(stringarray){
  single <- NULL
  for(i in seq(length(stringarray))){
    single <- paste(single,as.character(stringarray[i]),'')
  }
  return(single)
}
# =============================================================
# tree: tree structure from DataImportAndTreeGrown
# dat: data for covariates from DataImportAndTreeGrown
# cate: categorical covariates from DataImportAndTreeGrown
# link_group: TRUE/FALSE. distinguide the ribbon, when TRUE, target ribbons are distinguished for trt=T/C, red when target nodes recieve treatment treatment.effect.size>cost, grey when no treatment
# =============================================================
SankeyNetworkPlot <- function(tree,dat,cate,link_group){
  rownames(tree) <- seq(0,dim(tree)[1]-1,1)
  Link <- data.frame(source=character(),target=character(),IDsource=numeric(),IDtarget=numeric(),stringsAsFactors = FALSE)
  i <- 1
  for (j in seq(2,dim(tree)[1],1)) {
    nodeid <- tree[j,'node']
    isleft <- (substr(tree[j,'node'],start = nchar(tree[j,'node']),stop = nchar(tree[j,'node']))==1)
    node.parent.id <- substr(tree[j,'node'],start = 1,stop = nchar(tree[j,'node'])-1)
    node.parent <- tree[tree$node==node.parent.id,]
    var.name <- as.character(node.parent[,'vname'])
    if(is.element(var.name,cate)){
      cut <- as.character(unlist(strsplit(as.character(node.parent[,'cut']),split=" ")))
      cut.index <- is.element(levels(dat[,var.name]),cut)
      cut.val.left <- levels(dat[,var.name])[cut.index];   cut.val.left <- CombineLevelsString(cut.val.left)
      cut.val.right <- levels(dat[,var.name])[!cut.index]; cut.val.right <- CombineLevelsString(cut.val.right)
      target <- ifelse(isleft==TRUE,paste(var.name,'=',cut.val.left,sep = ''),paste(var.name,'=',cut.val.right,sep = ''))
      source <- ifelse(node.parent.id==1,'all',GetNodeName(tree,dat,cate,node.parent.id))
    } else{
      cut.val <- as.character(node.parent[,'cut'])
      target <- ifelse(isleft==TRUE,paste(var.name,'<=',cut.val,sep = ''),paste(var.name,'>',cut.val,sep = ''))
      source <- ifelse(node.parent.id==1,'all',GetNodeName(tree,dat,cate,node.parent.id))
    }
    Link[i,'source'] <- source; Link[i,'target'] <- target; Link[i,'value'] <- tree[j,'size']
    Link[i,'IDsource'] <- rownames(node.parent); Link[i,'IDtarget'] <- rownames(tree[j,])
    if(link_group==TRUE){
      Link[i,'group'] <- ifelse(tree[j,'trt']==1,'1','0')
    }
    i <- i+1
  }


  n_occur <- data.frame(table(Link$target))
  n_occur <- n_occur[n_occur$Freq > 1,]
  if(dim(n_occur)[1]!=0){
    for (i in seq(dim(n_occur)[1])){
      target <- n_occur[i,'Var1']
      target.data <- Link[Link$target==target,]
      k <- NULL
      for (j in seq(dim(target.data)[1])) {
        k <- paste(k,' ',sep = '')
        id.row <- rownames(target.data[j,])
        Link[id.row,'target'] <- paste(Link[id.row,'target'],k,sep = '')
        id.target <- Link[id.row,'IDtarget']
        if(is.element(id.target,Link$IDsource)){
          Link[Link$IDsource==id.target,'source'] <- Link[id.row,'target']
        }
      }
    }
  }


  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    name=c(as.character(Link$source),
           as.character(Link$target)) %>% unique()
  )


  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  Link$IDsource <- match(Link$source, nodes$name)-1
  Link$IDtarget <- match(Link$target, nodes$name)-1

  if(link_group==TRUE){
    nodes$group <- c('nodes.group')
    my_color <- 'd3.scaleOrdinal() .domain(["1", "0","nodes.group"]) .range(["red", "silver","silver"])'

    p <- sankeyNetwork(Links = Link, Nodes = nodes,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name",LinkGroup = 'group',
                       NodeGroup = 'group',colourScale = my_color,fontSize = 15,
                       sinksRight=FALSE)
    p
  } else{
    p <- sankeyNetwork(Links = Link, Nodes = nodes,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name",fontSize = 15,
                       sinksRight=FALSE)
    p
  }
}

# =================================================
# THE grow.INT() FUNCTION CONSTRUCTS A LARGE TREE
# =================================================
grow.INT <- function(dat, test=NULL,
                     split.var,ctg=NULL,
                     min.ndsz=20, pt=pt, max.depth=15,alpha,
                     mtry=length(split.var))
{
  out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- NULL
  list.nd <- list(dat);
  if (!is.null(test)) list.test <- list(test)
  name <- 1
  #i <- 1
  while (length(list.nd)!=0) {
    for (i in 1:length(list.nd)){
      print(i)
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){
        test0 <- NULL
        if (!is.null(test)) test0 <- list.test[[i]]
        split <- partition.INT(list.nd[[i]], test0, name[i], min.ndsz=min.ndsz,
                               pt=pt, split.var=split.var, ctg=ctg, max.depth=max.depth, alpha=alpha, mtry=mtry)
        out <- rbind(out, split$info)
        if (!is.null(split$left) && is.null(test)) {
          temp.list <- temp.test <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)

        } else if (!is.null(split$left) && !is.null(test) && !is.null(split$left.test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          temp.test <- c(temp.test, list(split$left.test, split$right.test))
        }
      }
    }
    list.nd <- temp.list; list.test <- temp.test; name <- temp.name
    temp.list <- temp.test <- temp.name <- NULL
  }
  out$node <- as.character(out$node)
  out <- out[order(out$node),]
  return(out)
}


grow.INT.ILL <- function(dat,
                         split.var,ctg=NULL,
                         min.ndsz, pt, max.depth,
                         mtry,
                         alpha,loc)
{
  out <- list.nd <- temp.list <- temp.name <- NULL
  list.nd <- list(dat);
  name <- 1
  #i <- 1
  while (length(list.nd)!=0) {
    for (i in 1:length(list.nd)){
      print(i)
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){
        split <- partition.INT.ILL(list.nd[[i]], name[i], min.ndsz=min.ndsz,
                                   pt=pt, split.var=split.var, ctg=ctg, max.depth=max.depth,mtry=mtry,alpha,loc)
        out <- rbind(out, split$info)
        if (!is.null(split$left)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
        }
      }
    }
    list.nd <- temp.list; name <- temp.name
    temp.list <-  temp.name <- NULL
  }
  out$node <- as.character(out$node)
  out <- out[order(out$node),]
  return(out)
}


grow.INT.ILL.OptimalDecision <- function(dat,
                                         split.var,ctg=NULL,
                                         min.ndsz, pt, max.depth,
                                         mtry,
                                         alpha,loc)
{
  out <- list.nd <- temp.list <- temp.name <- NULL
  list.nd <- list(dat);
  name <- 1
  #i <- 1
  while (length(list.nd)!=0) {
    for (i in 1:length(list.nd)){
      print(i)
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){
        split <- partition.INT.ILL.OptimalDecision(list.nd[[i]], name[i], min.ndsz=min.ndsz,
                                                   pt=pt, split.var=split.var, ctg=ctg,
                                                   max.depth=max.depth,mtry=mtry,alpha,loc)
        out <- rbind(out, split$info)
        if (!is.null(split$left)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
        }
      }
    }
    list.nd <- temp.list; name <- temp.name
    temp.list <-  temp.name <- NULL
  }
  out$node <- as.character(out$node)
  out <- out[order(out$node),]
  return(out)
}

# =================================================
# THE PARTITION INTINIAL LARGE TREE
# =================================================
partition.INT <- function(dat, test=NULL, name="1", min.ndsz=20, pt=pt,
                          split.var, ctg=NULL,
                          max.depth=15, alpha = 0.05, mtry=length(split.var))
{
  # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL VARIABLES.
  #browser()
  call <- match.call(); out <- match.call(expand = F)
  out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
  name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
  n <- nrow(dat);
  if (!is.null(test)) {n.test <- NA; score.test <- NA;}  ########## TEST SAMPLE ########
  var <- vname <- NA; cut <- NA; max.score <- threhold.score <- qt(1-alpha/2,df=n-2);
  trt <- dat$trt; y <- dat$y; vnames <- colnames(dat)
  # COMPUTE THE TREATMENT EFFECT IN CURRENT NODE
  #browser()
  trt.effect <- NA; n.1 <- sum(trt==1); n.0 <- n-n.1 ; n0 <- n*pt
  if (min(n.1, n.0) >= n0) {trt.effect <- mean(y)}
  # CONTROL THE MAX TREE DEPTH
  depth <- nchar(name)
  if (depth <= max.depth && n >= min.ndsz && min(n.1, n.0) >= n0) {
    m.try <- ifelse(is.null(mtry), length(split.var), mtry)
    for(i in sample(split.var, size=m.try, replace=F)) {
      x <- dat[,i]; v.name <- vnames[i]; temp <- sort(unique(x));
      if(length(temp) > 1) {
        if (is.element(i,ctg)) zcut <- power.set(temp)                                                  ############################ CLASS VARIABLE
        else zcut <- temp[-length(temp)]
        # print(i); print(temp); print(zcut)
        for(j in zcut) {
          score <- NA
          if (is.element(i,ctg)) {grp <- sign(is.element(x, j)); cut1 <- paste(j, collapse=" ")}      ############################ CLASS VARIABLE
          else  {grp <- sign(x <= j); cut1 <- as.character(j)}
          #browser()
          score <- ttest(dat, z=grp, min.ndsz, pt)
          # print(cbind(var=i, cut=j, score=score))
          if (!is.na(score) && score > max.score) {max.score <- score; var <- i; vname <- v.name; cut <- cut1; best.cut<-j}
        }}}}
  if (!is.null(test)) {
    n.test <- nrow(test); score.test <- NA;
    if (!is.na(var)) {
      if (is.element(var,ctg)) grp.test <- sign(is.element(test[,var], best.cut))                              ############################
      else  grp.test <- sign(test[,var] <= best.cut)
      score.test <- ttest(test, z=grp.test, n0=(n0/2))
      if (!is.na(score.test)){
        out$name.l <- name.l; out$name.r <- name.r
        out$left.test <- test[grp.test==1,  ]
        out$right.test <- test[grp.test==0,  ]
        if (is.element(var,ctg)) {                                                                               ############################
          out$left  <- dat[is.element(dat[,var], best.cut),]
          out$right <- dat[!is.element(dat[,var], best.cut), ]}
        else {
          out$left  <- dat[dat[,var]<= best.cut,]
          out$right <- dat[dat[,var]> best.cut, ]
        }
      }
      else {var <- NA; vname <- NA; cut <- NA;  max.score <- NA}
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,
                             var = var, vname=vname, cut= cut, score=ifelse(max.score==threhold.score, NA, max.score),
                             score.test, size.test=n.test)
    }
    else {
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,
                             var = NA, vname=NA, cut= NA, score=NA,
                             score.test=NA, size.test=n.test)
    }
  }	else {
    if (!is.na(var)) {
      out$name.l <- name.l; out$name.r <- name.r
      if (is.element(var,ctg)) {                                                                               ############################
        out$left  <- dat[is.element(dat[,var], best.cut),]
        out$right <- dat[!is.element(dat[,var], best.cut), ]}
      else {
        out$left  <- dat[dat[,var]<= best.cut,]
        out$right <- dat[dat[,var]> best.cut, ]
      }
      out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=dim(out$left)[1],n.r=dim(out$right)[1],
                             trt.effect=trt.effect, lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]], var = var, vname=vname, cut= cut,
                             score=ifelse(max.score==threhold.score, NA, max.score))
    }
    else{
      out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=NA, n.r=NA, trt.effect=trt.effect,
                             lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                             var=NA, vname=NA, cut=NA, score=NA)
    }
  }
  out
}



# =================================================
# THE PARTITION INTINIAL LARGE TREE to predict treatment effect,
# use increment of log-loglikelihood as the score of splits
# =================================================
partition.INT.ILL <- function(dat, name="1", min.ndsz, pt,
                              split.var, ctg=NULL,
                              max.depth, mtry=length(split.var),alpha,loc)
{
  # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL VARIABLES.
  #browser()
  call <- match.call(); out <- match.call(expand = F)
  out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
  name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
  n <- nrow(dat);
  var <- vname <- NA; cut <- NA; max.score <- 0;
  trt <- dat$trt; y <- dat$y; vnames <- colnames(dat)
  # COMPUTE THE TREATMENT EFFECT IN CURRENT NODE
  #browser()
  trt.effect <- ifelse(loc,sum(dat$weight*y)/sum(dat$weight),mean(y))
  #compute variance of this group treatment effect.
  mu  <- mean(y)
  predvar <- dat$predictive_var
  sigma <- (sum(predvar) + sum((y-mu)^2))/n

  n.1 <- sum(trt==1); n.0 <- n-n.1 ; n0 <- n*pt
  # CONTROL THE MAX TREE DEPTH
  depth <- nchar(name)
  if (depth <= max.depth && n >= min.ndsz && min(n.1, n.0) >= n0) {
    print(c(n,min.ndsz))
    m.try <- ifelse(is.null(mtry), length(split.var), mtry)
    for(i in sample(split.var, size=m.try, replace=F)) {
      x <- dat[,i]; v.name <- vnames[i]; temp <- sort(unique(x));
      if(length(temp) > 1) {
        if (is.element(i,ctg)) zcut <- power.set(temp)                                                  ############################ CLASS VARIABLE
        else zcut <- temp[-length(temp)]
        # print(i); print(temp); print(zcut)
        for(j in zcut) {
          score <- NA
          if (is.element(i,ctg)) {grp <- sign(is.element(x, j)); cut1 <- paste(j, collapse=" ")}      ############################ CLASS VARIABLE
          else  {grp <- sign(x <= j); cut1 <- as.character(j)}
          #browser()
          if (loc) { score <- WeightedIncrementLL(dat, grp, min.ndsz ,pt) - alpha}
          else { score <- IncrementLL(dat, grp, min.ndsz ,pt) - alpha }
          #print(cbind(var=i, cut=j, score=score,var=variance))
          if (!is.na(score) && score > max.score) {max.score <- score; var <- i; vname <- v.name; cut <- cut1; best.cut<-j}
        }}}}
  if (!is.na(var)) {
    out$name.l <- name.l; out$name.r <- name.r
    if (is.element(var,ctg)) {                                                                               ############################
      out$left  <- dat[is.element(dat[,var], best.cut),]
      out$right <- dat[!is.element(dat[,var], best.cut), ]}
    else {
      out$left  <- dat[dat[,var]<= best.cut,]
      out$right <- dat[dat[,var]> best.cut, ]
    }
    out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=dim(out$left)[1],n.r=dim(out$right)[1],
                           trt.effect=trt.effect, lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]], var = var, vname=vname, cut= cut,
                           score=max.score, se=sqrt(sigma)/sqrt(n))
    #print(variance)
  }
  else{
    out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=NA, n.r=NA, trt.effect=trt.effect,
                           lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                           var=NA, vname=NA, cut=NA, score=NA, se=sqrt(sigma)/sqrt(n))
    #print(variance)
  }
  return(out)
}


# =================================================
# THE PARTITION INTINIAL LARGE TREE to derive optimal decision z*
# use mean of prediction in each node
# =================================================
partition.INT.ILL.OptimalDecision <- function(dat, name="1", min.ndsz, pt,
                                              split.var, ctg=NULL,
                                              max.depth, mtry=length(split.var),alpha,loc)
{
  # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL VARIABLES.
  #browser()
  call <- match.call(); out <- match.call(expand = F)
  out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
  name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
  n <- nrow(dat);
  var <- vname <- NA; cut <- NA; max.score <- 0;
  trt <- dat$trt; y <- dat$optimal.prob; vnames <- colnames(dat);
  tau <- dat$y
  trt.effect <- ifelse(loc,sum(dat$weight*tau)/sum(dat$weight),mean(tau))

  #browser()
  mu  <- mean(y)

  n.1 <- sum(trt==1); n.0 <- n-n.1 ; n0 <- n*pt
  # CONTROL THE MAX TREE DEPTH
  depth <- nchar(name)
  if (depth <= max.depth && n >= min.ndsz && min(n.1, n.0) >= n0) {
    print(c(n,min.ndsz))
    m.try <- ifelse(is.null(mtry), length(split.var), mtry)
    for(i in sample(split.var, size=m.try, replace=F)) {
      x <- dat[,i]; v.name <- vnames[i]; temp <- sort(unique(x));
      if(length(temp) > 1) {
        if (is.element(i,ctg)) zcut <- power.set(temp)                                                  ############################ CLASS VARIABLE
        else zcut <- temp[-length(temp)]
        # print(i); print(temp); print(zcut)
        for(j in zcut) {
          score <- NA
          if (is.element(i,ctg)) {grp <- sign(is.element(x, j)); cut1 <- paste(j, collapse=" ")}      ############################ CLASS VARIABLE
          else  {grp <- sign(x <= j); cut1 <- as.character(j)}
          #browser()
          if (loc) { score <- WeightedIncrementLL(dat, grp, min.ndsz ,pt) - alpha}
          else { score <- IncrementLLOptimalDecision(dat, grp, min.ndsz) - alpha }
          #print(cbind(var=i, cut=j, score=score,var=variance))
          if (!is.na(score) && score > max.score) {max.score <- score; var <- i; vname <- v.name; cut <- cut1; best.cut<-j}
        }}}}
  if (!is.na(var)) {
    out$name.l <- name.l; out$name.r <- name.r
    if (is.element(var,ctg)) {                                                                               ############################
      out$left  <- dat[is.element(dat[,var], best.cut),]
      out$right <- dat[!is.element(dat[,var], best.cut), ]}
    else {
      out$left  <- dat[dat[,var]<= best.cut,]
      out$right <- dat[dat[,var]> best.cut, ]
    }
    out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=dim(out$left)[1],n.r=dim(out$right)[1],trt.effect=trt.effect,
                           prob.optimal=mu, lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]], var = var, vname=vname, cut= cut,
                           score=max.score)
    #print(variance)
  }
  else{
    out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=NA, n.r=NA, trt.effect=trt.effect, prob.optimal=mu,
                           lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                           var=NA, vname=NA, cut=NA, score=NA)
    #print(variance)
  }
  return(out)
}


# ==============================================
# PREDICT TREATMENT EFFECT
# tree: tree structure
# test: data for prediction
# ctg: categorical covariate index
# ==============================================
predict_initree_ITE <- function(tree,test,ctg){
  #browser()
  y <- NULL
  for (i in seq(dim(test)[1])){
    obs <- test[i,]
    node.id <- 1
    node <- tree[tree$node==node.id,]
    while (!is.na(node$vname)){
      varname <- as.character(node$vname)
      if(is.element(node$var,ctg)){
        cut.value <- as.character(node$cut)
        if(is.element(obs[,varname],strsplit(cut.value,' '))){
          node.id <- paste(node.id,'1',sep = '')
        } else{
          node.id <- paste(node.id,'2',sep = '')
        }
      } else{
        # cant use as.numeric(node$cut); below is efficient thant as.numeric(as.character(node$cut))
        cut.value <- as.numeric(levels(node$cut)[node$cut])
        if(obs[,varname]<=cut.value){
          node.id <- paste(node.id,'1',sep = '')
        }else{
          node.id <- paste(node.id,'2',sep = '')
        }
      }
      node <- tree[tree$node==node.id,]
    }
    y <- c(y,node$trt.effect)
  }
  return(y)
}


# ==============================================
# PREDICT PROBABILITY OF TREATMENT
# tree: tree structure
# test: data for prediction
# ctg: categorical covariate index
# ==============================================
predict_initree_OptimalDecision <- function(tree,test,ctg){
  #browser()
  y <- NULL
  for (i in seq(dim(test)[1])){
    obs <- test[i,]
    node.id <- 1
    node <- tree[tree$node==node.id,]
    while (!is.na(node$vname)){
      varname <- as.character(node$vname)
      if(is.element(node$var,ctg)){
        cut.value <- as.character(node$cut)
        if(is.element(obs[,varname],strsplit(cut.value,' '))){
          node.id <- paste(node.id,'1',sep = '')
        } else{
          node.id <- paste(node.id,'2',sep = '')
        }
      } else{
        # cant use as.numeric(node$cut); below is efficient thant as.numeric(as.character(node$cut))
        cut.value <- as.numeric(levels(node$cut)[node$cut])
        if(obs[,varname]<=cut.value){
          node.id <- paste(node.id,'1',sep = '')
        }else{
          node.id <- paste(node.id,'2',sep = '')
        }
      }
      node <- tree[tree$node==node.id,]
    }
    y <- c(y,node$prob.optimal)
  }
  return(y)
}


# ==============================================
# PREDICT node ID
# tree: tree structure
# test: data for prediction
# ctg: categorical covariate index
# ==============================================
predict_initree_NodeID <- function(tree,test,ctg){
  y <- c()
  for (i in seq(dim(test)[1])){
    obs <- test[i,]
    node.id <- 1
    node <- tree[tree$node==node.id,]
    while (!is.na(node$vname)){
      varname <- as.character(node$vname)
      if(is.element(node$var,ctg)){
        cut.value <- as.character(node$cut)
        if(is.element(obs[,varname],strsplit(cut.value,' '))){
          node.id <- paste(node.id,'1',sep = '')
        } else{
          node.id <- paste(node.id,'2',sep = '')
        }
      } else{
        # cant use as.numeric(node$cut); below is efficient thant as.numeric(as.character(node$cut))
        cut.value <- as.numeric(levels(node$cut)[node$cut])
        if(obs[,varname]<=cut.value){
          node.id <- paste(node.id,'1',sep = '')
        }else{
          node.id <- paste(node.id,'2',sep = '')
        }
      }
      node <- tree[tree$node==node.id,]
    }
    y <- c(y,node[,1])
  }
  return(y)
}

# ==============================================
# FUNCTION rdat() SIMULATES A SIMPLE DATA SET
# ==============================================

rdat <- function(n=100, K =50,
                 beta0=2, beta1=2, beta2=2, beta3=2, beta4=2, beta5=2,
                 sigma=1, cut1=.5, cut2=.5)
{
  trt <- sample(0:1, n, replace=T)
  #### Generate Covariates
  for (j in 1:4) assign(paste("x", j, sep=""),  sample(1:K, n, replace=T)/K)
  ###
  mean <- beta0 + beta1*trt + beta2*sign(x1<=cut1) + beta3*sign(x2<=cut2) +
    beta4*sign(x1<=cut1)*trt + beta5*sign(x2<=cut2)*trt
  y <- mean + rnorm(n, mean=0, sd=sigma)
  ##### Output
  data.frame(x1=x1, x2=x2, x3=x3, x4=x4, y=y, trt=trt)
}


# ===============================================================================
# FUNCTION ttest() COMPUTES THE t TEST FOR INTERACTION WITH CONTINUOUS RESPONSE
# ===============================================================================
ttest <- function(dat, z, min.ndsz, pt)
{
  #browser()
  trt <- dat$trt
  n.l.t <- sum((z==1)&(trt==1))/sum(z==1)
  n.l.c <- sum((z==1)&(trt==0))/sum(z==1)
  n.r.t <- sum((z==0)&(trt==1))/sum(z==0)
  n.r.c <- sum((z==0)&(trt==0))/sum(z==0)

  n <- nrow(dat)
  n.l <- sum(z==1)
  n.r <- n - n.l
  y.l <- dat[z==1,'y']
  y.r <- dat[z==0,'y']
  y <- dat$y
  score <- NA
  if (length(y)!=length(z)) stop("the vector z must have the same length as data.")
  t <- NA
  if ((min(n.l.t,n.l.c,n.r.t,n.r.c) >= pt)&&(n>min.ndsz)) {
    ybar.l <- mean(y.l); ybar.r <- mean(y.r)
    dif <- ybar.l - ybar.r
    t <- dif/sqrt(sum((y.l-ybar.l)^2)/((n.l-1)*n.l)+sum((y.r-ybar.r)^2)/((n.r-1)*n.r))
    score <- abs(t)
  }
  score
}


# ===============================================================================
# FUNCTION IncrementLL() COMPUTES Delta J(R_{m})
# ===============================================================================
IncrementLL <- function(dat, z, min.ndsz, pt)
{
  trt <- dat$trt
  n <- nrow(dat)
  n.l <- sum(z==1)
  n.r <- sum(z==0)
  # n.l.t <- sum((z==1)&(trt==1))/sum((z==1))
  # n.l.c <- sum((z==1)&(trt==0))/sum((z==1))
  # n.r.t <- sum((z==0)&(trt==1))/sum((z==0))
  # n.r.c <- sum((z==0)&(trt==0))/sum((z==0))
  #if ((n.l<min.ndsz)|(n.r<min.ndsz)|(min(n.l.t,n.l.c,n.r.t,n.r.c) < pt)) {
  if ((n.l<min.ndsz)|(n.r<min.ndsz)) {
    score <- NA
    return(score)
  }

  y.l <- dat[z==1,'y']
  y.r <- dat[z==0,'y']
  y <- dat$y

  mu.l <- mean(y.l)
  mu.r <- mean(y.r)
  mu   <- mean(y)

  predvar.l <- dat[z==1,'predictive_var']
  predvar.r <- dat[z==0,'predictive_var']
  predvar   <- dat$predictive_var

  sigma.l <- (sum(predvar.l) + sum((y.l-mu.l)^2))/n.l
  sigma.r <- (sum(predvar.r) + sum((y.r-mu.r)^2))/n.r
  sigma   <- (sum(predvar) + sum((y-mu)^2))/n

  ll.l <- -1/2*n.l*log(sigma.l)
  ll.r <- -1/2*n.r*log(sigma.r)
  ll   <- -1/2*n*log(sigma)

  score <- ll.l+ll.r-ll

  return(score)
}


WeightedIncrementLL <- function(dat, z, min.ndsz, pt)
{
  trt <- dat$trt
  # n <- nrow(dat)
  # n.l <- sum(z==1)
  # n.r <- sum(z==0)

  n <- sum(dat$weight)
  n.l <- sum(dat[z==1,'weight'])
  n.r <- sum(dat[z==0,'weight'])

  if ((n.l<min.ndsz)|(n.r<min.ndsz)) {
    score <- NA
    return(score)
  }

  y.l <- dat[z==1,'y']*dat[z==1,'weight']
  y.r <- dat[z==0,'y']*dat[z==0,'weight']
  y <- dat$y*dat$weight

  mu.l <- sum(y.l)/n.l
  mu.r <- sum(y.r)/n.r
  mu   <- sum(y)/n

  predvar.l <- dat[z==1,'predictive_var']*(dat[z==1,'weight'])^2
  predvar.r <- dat[z==0,'predictive_var']*(dat[z==0,'weight'])^2
  predvar   <- dat$predictive_var*(dat$weight)^2

  sigma.l <- (sum(predvar.l) + sum((y.l-mu.l)^2))/n.l
  sigma.r <- (sum(predvar.r) + sum((y.r-mu.r)^2))/n.r
  sigma   <- (sum(predvar) + sum((y-mu)^2))/n

  ll.l <- -1/2*n.l*log(sigma.l)
  ll.r <- -1/2*n.r*log(sigma.r)
  ll   <- -1/2*n*log(sigma)

  score <- ll.l+ll.r-ll

  return(score)
}


# ===============================================================================
# COMPUTES Delta J(R_{m}) to predict the optimal decision z*
# ===============================================================================
IncrementLLOptimalDecision <- function(dat, z, min.ndsz)
{
  n.l <- sum(z==1)
  n.r <- sum(z==0)

  if ((n.l<min.ndsz)|(n.r<min.ndsz)) {
    score <- NA
    return(score)
  }

  y.l <- dat[z==1,'optimal.prob']
  y.r <- dat[z==0,'optimal.prob']
  y <- dat[,'optimal.prob']

  mu.l <- mean(y.l)
  mu.r <- mean(y.r)
  mu   <- mean(y)

  scoreParent <- sum(y*log(mu)+(1-y)*log(1-mu))
  scoreLeft <- sum(y.l*log(mu.l)+(1-y.l)*log(1-mu.l))
  scoreRight <- sum(y.r*log(mu.r)+(1-y.r)*log(1-mu.r))

  score <- scoreLeft + scoreRight - scoreParent

  return(score)
}



# =================================
# METHOD I: THE TEST SAMPLE METHOD
# =================================

# -----------------------------------------------------------------------------
# THE PRUNING ALGORITHM GOOD FOR THE TREE SIZE SELECTION VIA TEST SAMPLE METHOD
# -----------------------------------------------------------------------------
#tre <- tree
prune.size.testsample <- function(tre)
{
  #browser()
  out <- as.list(NULL)
  ntest <- as.numeric(tre[1, ncol(tre)])
  if(is.null(dim(tre))) stop("No Need to Prune Further.")
  result <- NULL; n.tmnl <- sum(is.na(tre$var)); subtree <- 1
  a <- cbind(Ga.2=20, Ga.3=30, Ga.4=40, Ga.BIC=log(ntest))
  max.Ga <- rep(-1e20, 4); size <- rep(0, 4); btree <-as.list(1:4)
  while (n.tmnl > 1 ) {
    # print(tre)
    internal <- tre$node[!is.na(tre$cut)]; l <- length(internal);
    r.value <- 1:l
    for(i in 1:l) {
      branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
      score <- as.numeric(as.vector(branch$score))
      r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
    }
    alpha <- min(r.value)
    nod.rm <- internal[r.value == alpha];
    if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
    G <- sum(as.numeric(as.vector(tre$score)), na.rm=T);
    Ga <- G - a*l
    for (k in 1:4){if (Ga[k] > max.Ga[k]) {max.Ga[k] <- Ga[k]; size[k] <- n.tmnl; btree[[k]] <- tre}}
    result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre),
                                  size.tmnl=nrow(tre)-l, alpha=alpha, G=G, Ga))
    tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
    tre[match(nod.rm, tre$node), c("var", "vname", "cut", "score")] <- NA
    n.tmnl <- sum(is.na(tre$cut))
    if (n.tmnl ==1) {for (k in 1:4){if (0 > max.Ga[k]) {max.Ga[k] <- 0; size[k] <- 1; btree[[k]] <- tre}}}
    subtree <- subtree + 1
  }
  # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
  result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre),
                                size.tmnl=1, alpha=9999, G=0, Ga=cbind(Ga.2=0, Ga.3=0, Ga.4=0, Ga.BIC=0)))
  result <- as.data.frame(result)
  out$result <- result; out$size <- size; out$btree <- btree
  out
}



# ==========================================================
# FUNCTION de() FINDS ALL THE DESCENDENTS OF NODE x IN tree
# ==========================================================
de <- function(x, tree)
{
  #browser()
  if(length(x) != 1) stop("The length of x in function de must be 1.")
  y <- tree$node;  de <- NA
  if(sum(match(x, y), na.rm = T) != 0) {
    temp <- 1:length(y)
    start <- match(x, y) + 1
    end <- length(y)
    if(start <= length(y) & nchar(y[start]) > nchar(x)) {
      temp1 <- temp[temp >= start & nchar(y) <= nchar(x)][1] - 1  #find the subtree's the least right leaf node, then the next node is another tree, therefore finding the first node whose index larger the root node and depth is the same as the root node, then one node before this node is the end of the subtree (the last node of the subtree)
      if(!is.na(temp1)) end <- temp1
      de <- y[start:end]
    }}
  de
}



# ===========================================================================
# THE power.set() FUNCTION PROVIDES THE POWER SET FOR A CATEGORICAL VARIABLE
# ===========================================================================
power.set <- function(x) {
  if(length(x) == 0) return(vector(mode(x), 0))
  x <- sort(unique(x)); n <- length(x); K <- NULL
  for(m in x) K <- rbind(cbind(K, FALSE), cbind(K, TRUE))
  out <- apply(K, 1, function(x, s) s[x], s = x)
  out <- out[-c(1, length(out))]
  l <- length(out); i <- 1
  out[!sapply(out, length)>=ceiling(n/2+.5)]
}


# THIS FUNCTION as.numeric.factor() CONVERTS FACTOR INTO NUMERIC
as.numeric.factor <- function(x){as.numeric(levels(x))[x]}





# ========================================================================
# FUNCTION obtain.btree() OBTAINS THE BEST SUBTREE WITH KNOW SIZE bsize=
# ========================================================================
obtain.btree <- function(tre, bsize=6)
{
  btre <- NULL
  if (bsize==1) { btre <- tre[1,]; btre[, c("var", "cut", "score", "score.test")] <- NA}
  else if (bsize <1) stop("THE BEST TREE SIZE bsize= MUST BE >=1!")
  else {
    n.tmnl <- sum(is.na(tre$cut)); indicator <- T
    if (bsize > n.tmnl) stop("THE BEST TREE SIZE bsize PROVIDED IS LARGER THAN THE FULL TREE THAT YOU HAVE GROWN.")
    while (n.tmnl >= bsize && indicator ==T) {
      # print(tre); print(cbind(n.tmnl, bsize))
      internal <- tre$node[!is.na(tre$cut)]; l <- length(internal);
      r.value <- 1:l
      for(i in 1:l) {
        branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
        score <- as.numeric(as.vector(branch$score))
        r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
      }
      alpha <- min(r.value)
      nod.rm <- internal[r.value == alpha];
      tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
      tre[match(nod.rm, tre$node), c("var", "vname", "cut", "score", "score.test")] <- NA
      n.tmnl <- sum(is.na(tre$cut))
      # print(cbind(n.tmnl, bsize))
      if (n.tmnl==bsize) {btre <- tre; print(btre); indicator <- F}
    }
  }
  if (is.null(btre)) print(paste("The optimally-pruned subtree sequence does not have a subtree of bsize = ", bsize, sep=""))
  return(btre)
}



# ===================================================================
# PLOTTING IT TREE STRUCTURE, MODIFIED FROM PETER CALHOUN'S CODES
# ===================================================================
plot.tree <- function(tree, textDepth=3, lines="rectangle"){
  depth<-max(nchar(tree[,1]))
  par(xaxs='i')
  par(mar=c(1,1,1,1))
  par(xpd=TRUE)
  plot(1, type="n", xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), axes=FALSE,xaxs="i",yaxs="i")
  nodes<-tree$node
  nodesBin<-gsub("1", "0", nodes)
  nodesBin<-gsub("2", "1", nodesBin)
  lastObs<-nchar(nodesBin)
  nodesBin<-substr(nodesBin,2,lastObs)
  var <- tree$vname
  cut <- as.character(tree$cut)
  size <- tree$size

  for(i in 1:length(nodesBin)){
    nChar<-nchar(nodesBin[i])
    if(!is.na(var[i])){
      if(lines=="rectangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar)/(depth+1)))

        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))

        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))

      } else if(lines=="triangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))

        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
      }

      if(nChar <= textDepth){
        text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)+1/(depth+20),
             bquote(.(as.character(var[i]))<=.(cut[i])),cex=1)
      }
    } else {
      if(nChar <= textDepth){
        text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),(depth-nChar)/(depth+1),
             paste("n=", size[i], sep=""),cex=1, offset=1)
      }
    }
  }
}


# ================================================================
# TREE AUTOMATICALLY GENERATION; DATA IS CONVERTED INTO PARTY CLASS
# i: node (like 1,11,112)
# data: tree
# no return
# ================================================================
tree_search <- function(i,data){
  #browser()
  if(!is.element(i,data$node)){
    print(i)
    if (substr(i,nchar(i),nchar(i))==1){
      ii <- paste('split=sp_',substr(i,1,nchar(i)-1),sep='')
      ss <<- gsub(ii,'',ss,fixed=TRUE)
      return(NULL)
    } else{
      return(NULL)
    }
  }
  print(i)
  n <<- n+1
  left.kid <- paste(i,1,sep = '')
  right.kid <- paste(i,2,sep = '')

  if (substr(i,nchar(i),nchar(i))==1){
    if((!is.element(left.kid,data$node))&(!is.element(right.kid,data$node))){
      ss <<- paste(ss,'partynode(',n,'L,',paste('split=sp_',i,sep = ''),'info=',round(data[data$node==i,'trt.effect'],3),',kids=list{',sep = '')
    }else{
      ss <<- paste(ss,'partynode(',n,'L,',paste('split=sp_',i,sep = ''),',kids=list{',sep = '')
    }
  } else{
    if((!is.element(left.kid,data$node))&(!is.element(right.kid,data$node))){
      ss <<- paste(ss,',partynode(',n,'L,',paste('split=sp_',i,sep = ''),'info=',round(data[data$node==i,'trt.effect'],3),',kids=list{',sep = '')
    }else{
      ss <<- paste(ss,',partynode(',n,'L,',paste('split=sp_',i,sep = ''),',kids=list{',sep = '')
    }
  }

  tree_search(paste(i,1,sep = ''),data)
  tree_search(paste(i,2,sep = ''),data)
  ss <<- paste(ss,'})',sep = '')
}


# ================================================================
# GENERATE SP_NODE <- PARTYSPLIT(VAR,BREAKS(INDEX))
# dat: data for growing tree
# tree: tree
# ================================================================
generate_split <- function(dat,tree){
  #browser()
  cate <- NULL
  varname <- as.character(unique(tree$vname[!is.na(tree$vname)]))
  for (var in varname){
    if(class(dat[,var])=='factor'){
      cate <- cbind(cate,var)
    }
  }

  tree_split <- tree[!is.na(tree$var),]
  split_party <- NULL
  for (i in 1:dim(tree_split)[1]){
    name <- as.character(tree_split[i,'vname'])
    if (is.element(name,cate)){
      n <- length(levels(dat[,name]))
      index.cut <- 'c('
      cut.val <- unlist(strsplit(as.character(tree_split[i,'cut']),split=" "))
      for (j in 1:n){
        if (is.element(levels(dat[,name])[j],cut.val)) {index.cut <- ifelse(j==n,paste(index.cut,'1L)',sep = ''),paste(index.cut,'1L,',sep = ''))}
        else {index.cut <- ifelse(j==n,paste(index.cut,'2L)',sep = ''),paste(index.cut,'2L,',sep = ''))}
      }
      split_party <- paste(split_party,'sp_',tree_split[i,'node'],'<-','partysplit(',tree_split[i,'var'],'L,','index=',index.cut,');',sep = '')
    } else{
      split_party <- paste(split_party,'sp_',tree_split[i,'node'],'<-','partysplit(',tree_split[i,'var'],'L,','breaks=',tree_split[i,'cut'],');',sep = '')
    }
  }
  return(split_party)
}


# =============================================================
# PRUNE THE TREE BASED ON THE COST
# i: node, like 1, 11, 112
# data: tree
# =============================================================
prune_search <- function(i,data){
  #browser()
  if(!is.element(i,data$node)){
    print(i)
    return(NULL)
  }
  print(i)

  left.subtree <- prune_search(paste(i,1,sep = ''),data)
  right.subtree <- prune_search(paste(i,2,sep = ''),data)
  tree <- rbind(data[data$node==i,],left.subtree,right.subtree)

  if((!is.null(left.subtree))&(!is.null(right.subtree))){
    n.l <- dim(left.subtree)[1]; n.r <- dim(right.subtree)[1]
    if(n.l==n.r){
      n.t <- sum(left.subtree$trt == right.subtree$trt)
      if (n.t==n.l){
        tree <- data[data$node==i,]
        tree[,c('n.l','n.r','var','vname','cut','score')] <- NA
        return(tree)
      } else{
        return(tree)
      }
    } else{
      return(tree)
    }
  } else{
    return(tree)
  }
}

# ====================================
# EXTRACT AND PRINT RULES
# ====================================
print_leaf_nodes <- function(i,path,data,dat,cate){
  #browser()
  if(!is.element(i,data$node)){
    return(NULL)
  }

  if (i!=1){
    parent.node <- substr(i,1,nchar(i)-1)
    vname <- as.character(data[data$node==parent.node,'vname'])
    if(substr(i,nchar(i),nchar(i))==1){
      if(is.element(vname,cate)){
        cut <- as.character(unlist(strsplit(as.character(data[data$node==parent.node,'cut']),split=" ")))
        cut.index <- is.element(levels(dat[,vname]),cut)
        cut.val <- levels(dat[,vname])[cut.index]
        cut <- NULL
        for(j in 1:length(cut.val)){
          cut <- paste(cut,cut.val[j],sep = ' ')
        }
        vname.con <- paste(vname,'in', cut, sep = ' ')
        if (parent.node==1){
          path <- paste(path,'{',vname.con,sep = ' ')
        } else{
          path <- paste(path,vname.con,sep = ' AND ')
        }
      } else{
        cut <- data[data$node==parent.node,'cut']
        vname.con <- paste(vname,'<=', cut, sep = ' ')
        path <- ifelse(parent.node==1,paste(path,'{',vname.con,sep = ' '),paste(path,vname.con,sep = ' AND '))
      }
    } else{
      if(is.element(vname,cate)){
        cut <- as.character(unlist(strsplit(as.character(data[data$node==parent.node,'cut']),split=" ")))
        cut.index <- !is.element(levels(dat[,vname]),cut)
        cut.val <- levels(dat[,vname])[cut.index]
        cut <- NULL
        for(j in 1:length(cut.val)){
          cut <- paste(cut,cut.val[j],sep = ' ')
        }
        vname.con <- paste(vname,'in', cut, sep = ' ')
        if (parent.node==1){
          path <- paste(path,'{',vname.con,sep = ' ')
        } else{
          path <- paste(path,vname.con,sep = ' AND ')
        }
      } else{
        cut <- data[data$node==parent.node,'cut']
        vname.con <- paste(vname,'>', cut, sep = ' ')
        path <- ifelse(parent.node==1,paste(path,'{',vname.con,sep = ' '),paste(path,vname.con,sep = ' AND '))
      }
    }
  }

  if((!is.element(paste(i,1,sep = ''),data$node))&(!is.element(paste(i,2,sep = ''),data$node))){

    path <- paste(path,'}',
                  'ate=',round(data[data$node==i,'trt.effect'],3),
                  '(',round(data[data$node==i,'lower'],3),',',round(data[data$node==i,'upper'],3),')',
                  'size=',data[data$node==i,'size'],
                  't%=', round(data[data$node==i,'n.t']/data[data$node==i,'size'],3)*100,
                  'se=', format(data[data$node==i,'se'], digits=3),
                  sep = ' ')
    print(path)
  }

  print_leaf_nodes(paste(i,1,sep = ''),path,data,dat,cate)
  print_leaf_nodes(paste(i,2,sep = ''),path,data,dat,cate)

}

# ====================================
# compute the number of leaf node vs. sum of increment of loglikelihood
# tree is the initial large tree whose alpha=0
# alpha is a given pernalty term,
# we aim to plot two curves: x=number of leaves, y = 1) loglikeloood of the tree; 2) tree score (LL-alpha*#leaeves)
# ====================================
LLAndTreescoreGivenAlpha <- function(dat,tree,alpha){
  list.nd <- temp.list <- NULL
  list.nd <- tree[1,]

  n <- nrow(dat)
  y <- dat$y; mu<- mean(y)
  predvar   <- dat$predictive_var
  sigma   <- (sum(predvar) + sum((y-mu)^2))/n
  num.leaf <- 1
  LLTree  <- -1/2*n*log(sigma)
  Scoretree <- LLTree - num.leaf*alpha
  leaf_LL_Score <- c(num.leaf,LLTree,Scoretree)

  while (!is.null(list.nd)){
    for (i in 1:dim(list.nd)[1]){
      node <- list.nd[i,]
      if (is.element(paste(node[,1],'1',sep = ''),tree$node) && is.element(paste(node[,1],'2',sep = ''), tree$node)) {
        left <- as.numeric(paste(node[,1],'1',sep = '')); right <- as.numeric(paste(node[,1],'2',sep = ''))
        temp.list <- rbind(temp.list, tree[tree$node==left,], tree[tree$node==right,])
        num.leaf <- (num.leaf-1) + 2
        LLTree <- LLTree + node[,'score']
        Scoretree <- LLTree - num.leaf*alpha
        leaf_LL_Score <- rbind(leaf_LL_Score,c(num.leaf,LLTree,Scoretree))
      }
    }
    list.nd <- temp.list;
    temp.list <- NULL
  }
  leaf_LL_Score <- as.data.frame(leaf_LL_Score)
  colnames(leaf_LL_Score) <- c('treesize','LL','score')
  return(leaf_LL_Score)
}

# score = c(rank(delta LL(namely, score in tree) in each internal node, increasing))
# for (i in score){
#   continue split, step i, check object function C(sum of leaf nodes sigma - alpha*leaf_nodes), if Ci < Ci-1, then stop growing,
#   compute sum of log-likelihood of all the data of test sample
#   repeat 5 times,
#}
# compute alpha from bottom to top
# ======================================================
#tre <- tree
Get_Alpha_TreeSize <- function(tre) {
  n.tmnl <- sum(is.na(tre$cut));
  Alpha <- NULL
  TreeSize <- NULL
  while (n.tmnl > 1) {
    internal <- tre$node[!is.na(tre$cut)]; l <- length(internal);
    g.value <- 1:l
    for(i in 1:l) {
      branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
      score <- branch$score
      g.value[i] <- sum(score, na.rm=T) / (sum(is.na(score))-1)
    }
    alpha <- min(g.value)
    Alpha <- c(Alpha,alpha)
    nod.rm <- internal[g.value == alpha];
    tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
    tre[match(nod.rm, tre$node), c("var", "vname", "cut", "score")] <- NA
    n.tmnl <- sum(is.na(tre$cut))
    treesize <- sum(is.na(tre$score))
    TreeSize <- c(TreeSize,treesize)
    #print(treesize)
    #print(tre)
  }
  return(list(Alpha,TreeSize))
}


GetTrainSample <- function(samdat,n_sam){
  a <- sample(samdat,n_sam)
  return(a)
}


test.Rsquare <- function(test, tree,ctg){
  Outcome <- predict_initree_ITE(tree,test,ctg)
  SSres <- sum((test$y-Outcome)^2+test$predictive_var)
  SStot <- sum((dat$y-mean(dat$y))^2+dat$predictive_var)
  Rsquare <- 1 - SSres/SStot
  return(Rsquare)
}

test.logsigma <- function(test, tree,ctg){
  Outcome <- predict_initree_ITE(tree,test,ctg)
  Risk.test <- sum(log((test$y-Outcome)^2+test$predictive_var))
  return(Risk.test)
}

# input: follow up year, stage= (2,3,NULL)
# optput: list (processed data, xvars for tree growing)
DataImportForTreeGrow <- function(year,stage=NULL){
  if(year==5){
    data <- DataImportGeneral(year=5)
    load('~/Documents/iknl/code/policy_evaluation/post_y1_y0_5fp_update.RData')
  } else if (year==2){
    data <- DataImportGeneral(year=2)
    load('~/Documents/iknl/code/policy_evaluation/post_y1_y0_2fp_update.RData')
  } else{
    data <- DataImportGeneral(year=20156)
    load('~/Documents/iknl/code/policy_evaluation/post_y1_y0_2fp_20156.RData')
  }

  length_test <- (dim(post_y1_y0)[2])/2
  ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
  ate.test.mean <- apply(ate.test,2,mean)
  p1_ <- post_y1_y0[,c(1:length_test)]
  p_1 <- post_y1_y0[,c((length_test+1):(2*length_test))]
  expected.y1 <- apply(p1_, 2, mean)
  expected.y0 <- apply(p_1, 2, mean)

  data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                           ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
  data$ate <- ate.test.mean
  data$predictive_var <- apply(ate.test,2,var)

  if(!is.null(stage)){
    data <- filter(data,Stage==stage)
    xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num','year_diagnosis',
                'sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
                "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
    )
  } else{
    xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num','year_diagnosis',
                'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
                "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
    )
  }

  dat <- data[,c(xvars,'guideline2018','guideline2014')]
  dat$trt <- data$combined_chemo
  dat$y <- data$ate
  dat$predictive_var <- data$predictive_var

  return(list(dat,xvars,p1_,p_1))
}



DataImportForTreeGrow20156 <- function(data,stage=NULL){
  data <- data_process_all_variable_update()
  data <- data_process_to_dummy_update(data)
  data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
  data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
  target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis>2014)
  data <- data[target_index,]

  load("post_y1_p0_2fp_20156.RData")
  length_test <- (dim(post_y1_y0)[2])/2
  ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
  ate.test.mean <- apply(ate.test,2,mean)

  data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                             ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
  data$ate <- ate.test.mean
  data$predictive_var <- apply(ate.test,2,var) #var is unbiased estimator of variance

  length_test <- (dim(post_y1_y0)[2])/2
  ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
  ate.test.mean <- apply(ate.test,2,mean)
  p1_ <- post_y1_y0[,c(1:length_test)]
  p_1 <- post_y1_y0[,c((length_test+1):(2*length_test))]
  expected.y1 <- apply(p1_, 2, mean)
  expected.y0 <- apply(p_1, 2, mean)



  if(!is.null(stage)){
    data <- filter(data,Stage==stage)
    xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num',
                'sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
                "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
    )
  } else{
    xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num',
                'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
                "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
    )
  }

  dat <- data[,c(xvars,'guideline2018','guideline2014')]
  dat$trt <- data$combined_chemo
  dat$y <- data$ate
  dat$predictive_var <- data$predictive_var

  return(list(dat,xvars,p1_,p_1))
}


OptimalDecisionGivenLoss <- function(rho,p1_,p_1,cost){
  rho <- rho
  p11 <- rho*sqrt(p1_*(1-p1_)*p_1*(1-p_1)) + p1_*p_1
  p10 <- p1_ - p11
  p01 <- p_1 - p11
  p00 <- 1- p11 - p10 - p01
  costControl <- cost[1,]
  costTreatment <- cost[2,]
  deltaCost <- costTreatment-costControl
  DeltaPosteriorCost <- deltaCost[1]*p00 + deltaCost[2]*p01 + deltaCost[3]*p10 + deltaCost[4]*p11
  OptimalDecisionProb <- apply(DeltaPosteriorCost<0, 2, mean)
  return(OptimalDecisionProb)
}


DataImport20156 <- function(){
  set.seed(123)
  xvars.all  <- c('age_at_diagnosis','lymph_assessed',
                  'male',#'female',
                  'BMI',
                  'MSI_0','MSI_1','MSI_unknown',#'MSI_missing',
                  'colon_perforation_0',"colon_perforation_1","colon_perforation_unknown",
                  #"colon_perforation_missing",
                  'lymphomatic_invation_no',
                  'lymphomatic_invation_yes',
                  'lymphomatic_invation_suspect',
                  'lymphomatic_invation_unknown',
                  #'lymphatic_vessel_invasion_missing'
                  'agio_invation_no',
                  'agio_invation_EMVI',
                  'agio_invation_IMVI',
                  'agio_invation_suspect',
                  #'agio_invation_NA',
                  'agio_invation_unknown',
                  "Grade_good_moderate","Grade_poor_undifferentiated",
                  #'pT0',
                  'pT1','pT2','pT3',#'pT4',
                  'pN0','pN1','pN2','pNX',#'pN.missing',
                  'cM0','cM1','cMX',#'cM.missing',
                  'max_LoS',
                  'Morphology_1','Morphology_2','Morphology_3','Morphology_4',
                  'Stage2','Stage3',
                  'highrisk_false','highrisk_true','highrisk_unknown', #Morphology_5
                  'proximal','distal',
                  'BRAF_0','BRAF_1','BRAF_9',
                  'RAS_0','RAS_1','RAS_9',
                  'ASAclass_1','ASAclass_2','ASAclass_3','ASAclass_4','ASAclass_5','ASAclass_9',
                  'PS_notregistered','PS_unknown','PS_WHO1','PS_WHO2','PS_WHO3','PS_WHO4'
  )
  if(file_test("-f", "post_y1_p0_2fp_20156.RData")){
    data <- data_process_all_variable_update()
    data <- data_process_to_dummy_update(data)
    data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
    data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis>2014)
    data <- data[target_index,]

    load("post_y1_p0_2fp_20156.RData")
    length_test <- (dim(post_y1_y0)[2])/2
    ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
    ate.test.mean <- apply(ate.test,2,mean)
    #target population adding ate for each observation
    data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                             ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
    data$ate <- ate.test.mean
    data$predictive_var <- apply(ate.test,2,var) #var is unbiased estimator of variance
  } else{
    ntree <- 50
    data <- data_process_all_variable_update()
    data <- data_process_to_dummy_update(data)
    data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
    data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis>2014)
    data <- data[target_index,]

    set.seed(123)
    length_test <- dim(data)[1]
    train <- data[,c(xvars.all)]
    y <- data$combined_chemo
    test <- data[,xvars.all]
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    ps <- apply(post$prob.test,2,mean)
    data$bartps <- ps

    year <- 2
    selection <- data
    selection <- filter(selection,((vitfup<year)&(vitstat==1))|(vitfup>=year))
    selection[selection$vitfup>=year,"vitstat"] <- 0

    length_test <- dim(data)[1]
    train <- selection[,c('combined_chemo',xvars.all,'bartps')]
    y <- selection$vitstat
    test <- data[,c('combined_chemo',xvars.all,'bartps')]
    test_bart_1 <- test_bart_0 <- test
    test_bart_1$combined_chemo <- 1
    test_bart_0$combined_chemo <- 0
    test <- rbind(test_bart_1,test_bart_0)
    post <- pbart(train,y,test,ntree = ntree)
    post$prob.test <- pnorm(post$yhat.test)
    post$prob.test <- 1- post$prob.test

    post_y1_y0 <- post$prob.test
    save(post_y1_y0, file = "post_y1_p0_2fp_20156.RData")
    length_test <- (dim(post_y1_y0)[2])/2
    ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
    ate.test.mean <- apply(ate.test,2,mean)


    #target population adding ate for each observation
    data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                             ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
    data$ate <- ate.test.mean
    data$predictive_var <- apply(ate.test,2,var)
  }
  # =============================================================
  # POLICY FOR INTERPRETABLE TREE
  # =============================================================
  xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num',
              'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
              "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
  )
  dat <- data[,xvars]
  dat$trt <- data$combined_chemo
  dat$y <- ate.test.mean
  dat$predictive_var <- data$predictive_var
  return(list(dat,xvars))
}



# DataImportLocal <- function(year,ID){
#   data <- data_process_all_variable_update()
#   data <- data_process_to_dummy_update(data)
#
#   if(year==5){
#     target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis<2014)
#     data <- data[c(target_index),]
#     load('post_y1_y0_5fp_update.RData')
#   } else{
#     target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))
#     data <- data[c(target_index),]
#     load('post_y1_y0_2fp_update.RData')
#   }
#
#
#   length_test <- (dim(post_y1_y0)[2])/2
#   ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
#
#   return(ate.test[,ID])
# }


PlotAte <- function(ate){
  ate <- as.data.frame(ate)
  colnames(ate) <- c("ITE")
  p <- ggplot(data=ate,aes(ITE))+
       geom_density(color="lightblue", fill="lightblue") +
       ggtitle("Distribution of x_{i} individual treatment effect")
  return(p)
}


# input:
# alpha: penalty term of mode complexity(number of internal splits)
# alpha=0: to grow a large inital tree without prunning
# alpha >0: to grow and prune a tree given a fixed alpha
# dat(processed dat),xvars(for tree grow),min.ndsz(minimun number of obs in node),
# pt(minimun overlap between treatment and control groups), max.depth(max depth of tree)
TreeGrownAndPrune <- function(dat,xvars,min.ndsz,pt,max.depth,alpha,loc){
  # COLUMNS OF COVARIATES
  split.var <- 1:length(xvars);ctg=6:length(xvars); mtry=length(split.var)
  # FLEXIBILTIY OF TREE STRUCTURE
  min.ndsz=min.ndsz # mininum number of observations in a node
  pt=pt # minimun number of treatment/control in a node
  max.depth=max.depth # max depth of tree

  tree <- grow.INT.ILL(dat,
                       split.var,ctg,
                       min.ndsz, pt, max.depth,
                       mtry=length(split.var),alpha,loc)
  return(tree)
}


#' @param min.ndsz stopping rule, the minimum number observations in leaf node
#' @param pt stopping rule, the minimum value of percentage of patients getting treatment, minimum propensity score
#' @param max.depth stopping rule, the maximum depth of the tree
TreeGrownAndPruneOptimalDecision <- function(dat,xvars,min.ndsz,pt,max.depth,alpha,loc){
  # COLUMNS OF COVARIATES
  split.var <- 1:length(xvars);ctg=6:length(xvars); mtry=length(split.var)
  # FLEXIBILTIY OF TREE STRUCTURE
  min.ndsz=min.ndsz # mininum number of observations in a node
  pt=pt # minimun number of treatment/control in a node
  max.depth=max.depth # max depth of tree

  tree <- grow.INT.ILL.OptimalDecision(dat,split.var,ctg,min.ndsz, pt,
                                       max.depth,mtry=length(split.var),alpha,loc)
  return(tree)
}


# input: processed data, tree
GetCategoryVarInTree <- function(dat,tree){
  varname <- as.character(unique(tree$vname[!is.na(tree$vname)]))
  cate <- NULL
  for (var in varname){
    if(class(dat[,var])=='factor'){
      cate <- cbind(cate,var)
    }
  }
  return(cate)
}

# input: processed data,xvars for growing tree, initial large tree
# output: alpha and error of pruned tree
TreeComplexity5foldCV <- function(dat,xvars,tree,min.ndsz,max.depth,pt,loc){
  ##### generate test index for k-fold cross validation
  k <- 5
  n_test <- floor(nrow(dat)/k)
  Test <- list()
  Test[[1]] <- GetTrainSample(rownames(dat),n_test)
  Test[[2]] <- GetTrainSample(rownames(dat[!is.element(rownames(dat),Test[[1]]),]),n_test)
  Test[[3]] <- GetTrainSample(rownames(dat[!is.element(rownames(dat),c(Test[[1]],Test[[2]])),]),n_test)
  Test[[4]] <- GetTrainSample(rownames(dat[!is.element(rownames(dat),c(Test[[1]],Test[[2]],Test[[3]])),]),n_test)
  Test[[5]] <- GetTrainSample(rownames(dat[!is.element(rownames(dat),c(Test[[1]],Test[[2]],Test[[3]],Test[[4]])),]),n_test)

  Alpha_TreeSize <- Get_Alpha_TreeSize(tree)
  Alpha <- Alpha_TreeSize[[1]]
  TreeSize <- Alpha_TreeSize[[2]]
  Risk <- NULL
  split.var <- 1:length(xvars);ctg=6:length(xvars); mtry=length(split.var)

  CV <- function(i){
    test <- dat[is.element(rownames(dat),Test[[i]]),]
    train <- dat[!is.element(rownames(dat),Test[[i]]),]
    R.fold <- NULL

    for (j in 1:length(Alpha)){
      tree <- grow.INT.ILL(train,
                           split.var,ctg,
                           min.ndsz=min.ndsz, pt=pt, max.depth,
                           mtry,
                           alpha = Alpha[j],loc)
      R.test <- test.logsigma(test,tree,ctg)
      R.fold <- c(R.fold,R.test)
    }
    return(R.fold)
  }

  Risk <- mclapply(seq(k),CV,mc.cores = 3)
  Risk <- t(as.data.frame(Risk))
  Risk <- apply(Risk, 2, sum)/dim(dat)[1]
  Risk <- rbind(Alpha,TreeSize,Risk)

  return(Risk)
}

plotLeafScore <- function(df){
  df <- as.data.frame(df)
  p1 <- ggplot(data=df, aes(x=treesize, y=score, color=Alpha)) +
    geom_line(linetype = "solid")+
    geom_point() +
    ggtitle(expression("Tree size versus G(T;"~alpha~")"))+
    xlab("Tree Size")+
    ylab(expression("G(T;"~alpha~")"))+
         #color=expression(alpha),size=6)+
    theme(legend.position = c(0.8, 0.3),
          legend.text = element_text(size=10),
          legend.background = element_rect(fill=alpha('white', 0.5)),
          plot.title = element_text(size=10)) +
    scale_colour_discrete(expression(alpha))
  return(p1)
}

TreeSize_Risk <- function(df,vadjust,x_break){
  df <- df[,order(df[2,])]
  df <- t(df)
  df <- as.data.frame(df)
  df[,1] <- round(df[,1],digits = 1)
  n <- dim(df)[1]
  alpha.yposition <- max(df[,'Risk'])
  # if((floor(n)/2)<8){
  #   df[c(1,seq(3,floor(n)/2,1),seq(floor(n)/2+2,n)),1] <- NA
  # } else if((floor(n)/3)<8){
  #   df[c(1,3,4,5,6,seq(8,floor(n)/2,1),seq(floor(n)/2+2,n-1)),1] <- NA
  # } else{
  #   df[c(1,3,4,5,6,seq(8,floor(n)/3,1),seq(floor(n)/3+2,2*floor(n)/3,1),seq(2*floor(n)/3+2,n-1,1)),1] <- NA
  # }

  df[c(seq(2,floor(n)/3,1),seq(floor(n)/3+2,2*floor(n)/3,1),seq(2*floor(n)/3+2,n-1,1)),1] <- NA

  p1 <- ggplot(data=df, aes(x=TreeSize, y=Risk,label=Alpha)) +
    geom_point()+
    scale_x_continuous(breaks=seq(min(df[,2]),max(df[,2]),x_break)) +
    geom_line(linetype = "solid")+
    coord_cartesian(clip = 'off') +
    labs(#title = expression("Cross validated estimate of risk"),
         subtitle= expression(alpha),
         x="Tree Size",y="J(f2)")+
    theme(plot.margin = unit(c(1.5,1,1,1), "lines"),
          plot.subtitle = element_text(hjust = 0.5,vjust = 4.5,size=10),
          plot.title = element_text(vjust = 5.5,size=12),
          text = element_text(size=15)) +
    geom_text(y=max(alpha.yposition)+vadjust,size=3)
  return(p1)
}


PlotCostSensTree <- function(tree.pruned,dat,cate,cost){
  tree.pruned$trt <- 1*(tree.pruned$trt.effect > cost)
  tree.pruned <- prune_search(1,tree.pruned)
  path <- NULL
  print_leaf_nodes(1,path,tree.pruned,dat,cate)
  SankeyNetworkPlot(tree.pruned,dat,cate,link_group = TRUE)
}


Permute_VarIP <- function(i){
  n_sam <- nrow(dat)*0.9
  train.id <- sample(1:nrow(dat),n_sam)
  train <- dat[train.id,]
  test <- dat[-train.id,]
  tree <- TreeGrownAndPrune(train,xvars,min.ndsz=min.ndsz,pt=pt,max.depth=max.depth,alpha=0,loc=FALSE)
  #2. remaining test sample to compute G(T)
  ctg <- 6:length(xvars)
  score.bench <- test.logsigma(test,tree,ctg)
  #3. permute x_{j} then compute Gj(T), and compute G(T)-Gj(T)
  importance.permute <- NULL
  for (xvar in xvars) {
    test.permute <- test
    permute.var <- permute(test[,xvar])
    test.permute[,xvar] <- permute.var
    score.permute <- test.logsigma(test.permute,tree,ctg)
    importance.var <- (score.bench - score.permute)/score.bench
    importance.permute <- c(importance.permute,importance.var)
  }
  return(importance.permute)
}

Permute_VarIP_LC <- function(i){
  n_sam <- nrow(neighbour)*0.9
  train.id <- sample(1:nrow(neighbour),n_sam)
  train <- neighbour[train.id,]
  test <- neighbour[-train.id,]
  tree <- TreeGrownAndPrune(neighbour,xvars,min.ndsz=min.ndsz,pt=0,max.depth=max.depth,alpha=0,loc=FALSE)
  #2. remaining test sample to compute G(T)
  ctg <- 6:length(xvars)
  score.bench <- test.test.logsigma(test,tree,ctg)
  #3. permute x_{j} then compute Gj(T), and compute G(T)-Gj(T)
  importance.permute <- NULL
  for (xvar in xvars) {
    test.permute <- test
    permute.var <- permute(test[,xvar])
    test.permute[,xvar] <- permute.var
    score.permute <- test.logsigma(test.permute,tree,ctg)
    importance.var <- (score.bench - score.permute)/score.bench
    importance.permute <- c(importance.permute,importance.var)
  }
  return(importance.permute)
}


Treesize_Score <- function(R.ave,dat,tree,alpha.costum){
  alpha.cv <- R.ave[1,which(R.ave[3,]==max(R.ave[3,]))]; alpha.cv <- round(alpha.cv,digits = 1)
  alpha.aic <- 2
  alpha.bic <- log(nrow(dat))
  leaf_LL_Score.BIC <- LLAndTreescoreGivenAlpha(dat,tree,alpha.bic)
  leaf_LL_Score.AIC <- LLAndTreescoreGivenAlpha(dat,tree,alpha.aic)
  leaf_LL_Score.CV  <-  LLAndTreescoreGivenAlpha(dat,tree,alpha.cv)
  leaf_LL_Score.CUS <-  LLAndTreescoreGivenAlpha(dat,tree,alpha.costum)

  leaf_LL_Score.BIC$Alpha <- "BIC:log(n) "
  leaf_LL_Score.AIC$Alpha <- "AIC:2"
  leaf_LL_Score.CV$Alpha  <- paste('CV:',alpha.cv,sep = '')
  leaf_LL_Score.CUS$Alpha <- round(alpha.costum,digits = 1)
  leaf_LL_Score <- rbind(leaf_LL_Score.BIC,leaf_LL_Score.AIC,leaf_LL_Score.CV,leaf_LL_Score.CUS)
  pGalpha <- plotLeafScore(leaf_LL_Score)
  return(pGalpha)
}




PlotFI <- function(bootstrap.IP,xvars){
  bootstrap.IP <- as.data.frame(bootstrap.IP)
  bootstrap.IP <- apply(bootstrap.IP,1, mean)
  bootstrap.IP <- as.data.frame(bootstrap.IP)
  bootstrap.IP<-  as.data.frame(cbind(xvars,bootstrap.IP))
  colnames(bootstrap.IP) <- c("var",'importance')
  # convert factor to numeric: factor->character->numeric
  bootstrap.IP$importance <- as.numeric(as.character(bootstrap.IP$importance))
  bootstrap.IP <- filter(bootstrap.IP,importance!=0)
  # order var based on importance ranking
  bootstrap.IP$var <- with(bootstrap.IP, factor(var,
                                                levels=var[order(importance)]))
  pFI <- ggplot(data=bootstrap.IP, aes(x=var, y=importance)) +
    geom_bar(stat="identity") + coord_flip()
    #ggtitle("Feature importance via random forest")

  return(pFI)
}

ChooseAlpha <- function(R.ave,threhold){
  n <- dim(R.ave)[2]
  Rsqure <- 0
  stop <- FALSE
  while((!stop)&(n>=1)){
    if(abs(R.ave[3,n]-Rsqure)>=threhold){
      Rsqure <- R.ave[3,n]
      alpha <- R.ave[1,n]
      n <- n-1
    } else{
      stop <- TRUE
    }
  }
  return(alpha)
}



PolicyEvaluationDataImport <- function(year, load_data){
  # =============================================================
  # PREPARATION: IMPORT DATA, COMPUTE GUIDELINE, SELECT TARGET POPULATION, AND POTENTIAL OUTCOME COMPUTATIONS
  # =============================================================
  #browser()
  data <- data_process_all_variable_update()
  data <- data_process_to_dummy_update(data)
  data$guideline2018 <- getPrescribedTreatment_guideline2018(data)
  data$guideline2014 <- getPrescribedTreatment_guideline2014(data)
  # target population index
  if(year==5){
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))&(data$year_diagnosis<2014)
    data <- data[target_index,]
  } else{
    target_index <- (!((data$cM==1)|(data$cM=='1A')|(data$cM=='1B')))&((data$Stage==2)|(data$Stage==3))
    data <- data[target_index,]
  }

  file.name <- ifelse(year==5,'post_y1_y0_5fp_update.RData','post_y1_y0_2fp_update.RData')
  load(file.name)

  length_test <- (dim(post_y1_y0)[2])/2
  ate.test <- post_y1_y0[,1:length_test]-post_y1_y0[,(length_test+1):(2*length_test)]
  ate.test.mean <- apply(ate.test,2,mean)
  p1_ <- post_y1_y0[,c(1:length_test)]
  p_1 <- post_y1_y0[,c((length_test+1):(2*length_test))]
  expected.y1 <- apply(p1_, 2, mean)
  expected.y0 <- apply(p_1, 2, mean)

  #combine lymphomatic invation and agio invation together
  data$lym_agioinvasie <- as.factor(ifelse((data$lymphomatic_invation_yes==1)|(data$agio_invation_EMVI==1)|(data$agio_invation_IMVI==1),'yes',
                                           ifelse((data$lymphomatic_invation_no==1)&(data$agio_invation_no==1),'no','other')))
  data$ate <- ate.test.mean
  data$y1 <- expected.y1
  data$y0 <- expected.y0

  return(data)
}


ComputeLossOutValuePolicies <- function(year,data,tree.ate.hat,p1_,p_1,rho,cd,ct_seq){
  filename <- ifelse(length(ct_seq)==1,paste('loss_outcome_value_rho',rho,'_cost',ct_seq,'_',year,'fp.RData',sep = ''),
                     paste('loss_outcome_value_rho',rho,'_',year,'fp.RData',sep=''))
  # =============================================================
  # correlation of Y(1),Y(0): rho
  # theta
  # =============================================================
  rho <- rho
  p11 <- rho*sqrt(p1_*(1-p1_)*p_1*(1-p_1)) + p1_*p_1
  p10 <- p1_ - p11
  p01 <- p_1 - p11
  p00 <- 1- p11 - p10 - p01

  expected.y1 <- apply(p1_, 2, mean); var.y1 <- apply(p1_,2,var)
  expected.y0 <- apply(p_1, 2, mean); var.y0 <- apply(p_1,2,var)

  L0 <- cd*p00+cd*p10; expected.L0 <- apply(L0, 2, mean); var.L0 <- apply(L0, 2, var)
  V0 <- p_1 - L0; expected.V0 <- apply(V0, 2, mean); var.V0 <- apply(V0, 2, var)
  # =============================================================
  # COMPUTE POSTERIOR EXPECTED LOSS AND POSTERIOR EXPECTED OUTCOME
  # =============================================================
  guideline2018.action <- ifelse(data$guideline2018=='unknown',ifelse(data$y>0,1,0),data$guideline2018)

  # =============================================================
  loss <- data.frame(cost=numeric(),estimator=character(),mean=numeric(),
                     lower=numeric(),upper=numeric(),standard_error=numeric(),
                     QL=numeric(),QU=numeric(),stringsAsFactors = FALSE)
  outcome <- data.frame(cost=numeric(),estimator=character(),mean=numeric(),
                        lower=numeric(),upper=numeric(),standard_error=numeric(),
                        QL=numeric(),QU=numeric(),stringsAsFactors = FALSE)
  value <- data.frame(cost=numeric(),estimator=character(),mean=numeric(),
                      lower=numeric(),upper=numeric(),standard_error=numeric(),
                      QL=numeric(),QU=numeric(),stringsAsFactors = FALSE)

  i <- 0
  for (ct in ct_seq){ # from 0 to 1 (include 0 an 1)
    threshold <- ct/cd
    L1 <- (ct+cd)*p00+(ct+cd)*p01+ct*p10+ct*p11; expected.L1 <- apply(L1, 2, mean); var.L1 <- apply(L1, 2, var)
    V1 <- p1_- L1; expected.V1 <- apply(V1, 2, mean);  var.V1 <- apply(V1, 2, var);

    chemo.outcome <- expected.y1
    chemo.loss <- expected.L1
    chemo.value <- expected.V1
    chemo.outcome.var <- var.y1
    chemo.loss.var <- var.L1
    chemo.value.var <- var.V1

    surgery.outcome <- expected.y0
    surgery.loss <- expected.L0
    surgery.value <- expected.V0
    surgery.outcome.var <- var.y0
    surgery.loss.var <- var.L0
    surgery.value.var <- var.V0

    guideline2018.outcome <- ifelse(guideline2018.action==1,expected.y1,expected.y0)
    guideline2018.loss <- ifelse(guideline2018.action==1,expected.L1,expected.L0)
    guideline2018.value <- ifelse(guideline2018.action==1,expected.V1,expected.V0)
    guideline2018.outcome.var <- ifelse(guideline2018.action==1,var.y1,var.y0)
    guideline2018.loss.var <- ifelse(guideline2018.action==1,var.L1,var.L0)
    guideline2018.value.var <- ifelse(guideline2018.action==1,var.V1,var.V0)

    optimal.action <- ifelse(data$y>threshold,1,0)
    #optimal.action.r <- t(replicate(1000,optimal.action))
    optimal.outcome <- ifelse(optimal.action==1,expected.y1,expected.y0)
    optimal.loss <- ifelse(optimal.action==1,expected.L1,expected.L0)
    optimal.value <- ifelse(optimal.action==1,expected.V1,expected.V0)
    optimal.outcome.var <- ifelse(optimal.action==1,var.y1,var.y0)
    optimal.loss.var <- ifelse(optimal.action==1,var.L1,var.L0)
    optimal.value.var <- ifelse(optimal.action==1,var.V1,var.V0)


    tree.action <- ifelse(tree.ate.hat>threshold,1,0)
    tree.outcome <- ifelse(tree.action==1,expected.y1,expected.y0)
    tree.loss <- ifelse(tree.action==1,expected.L1,expected.L0)
    tree.value <- ifelse(tree.action==1,expected.V1,expected.V0)
    tree.outcome.var <- ifelse(tree.action==1,var.y1,var.y0)
    tree.loss.var <- ifelse(tree.action==1,var.L1,var.L0)
    tree.value.var <- ifelse(tree.action==1,var.V1,var.V0)


    i <- i+1
    loss[i,] <- append.data(chemo.loss,'chemo',ct,chemo.loss.var)
    outcome[i,] <- append.data(chemo.outcome,'chemo',ct,chemo.outcome.var)
    value[i,] <- append.data(chemo.value,'chemo',ct,chemo.value.var)

    i <- i+1
    loss[i,] <- append.data(surgery.loss,'surgery',ct,surgery.loss.var)
    outcome[i,] <- append.data(surgery.outcome,'surgery',ct,surgery.outcome.var)
    value[i,] <- append.data(surgery.value,'surgery',ct,surgery.value.var)

    i <- i+1
    loss[i,] <- append.data(guideline2018.loss,'guideline2018',ct,guideline2018.loss.var)
    outcome[i,] <- append.data(guideline2018.outcome,'guideline2018',ct,guideline2018.outcome.var)
    value[i,] <- append.data(guideline2018.value,'guideline2018',ct,guideline2018.value.var)

    i <- i+1
    loss[i,] <- append.data(optimal.loss,'optimal',ct,optimal.loss.var)
    outcome[i,] <- append.data(optimal.outcome,'optimal',ct,optimal.outcome.var)
    value[i,] <- append.data(optimal.value,'optimal',ct,optimal.value.var)

    i <- i+1
    loss[i,] <- append.data(tree.loss,'tree',ct,tree.loss.var)
    outcome[i,] <- append.data(tree.outcome,'tree',ct,tree.outcome.var)
    value[i,] <- append.data(tree.value,'tree',ct,tree.value.var)

    print(ct)
  }

  if(length(ct_seq)==1){
    loss <- cbind(optimal.loss,tree.loss,guideline2014.loss,guideline2018.loss,chemo.loss,surgery.loss,observation.loss)
    outcome <- cbind(optimal.outcome,tree.outcome,guideline2014.outcome,guideline2018.outcome,chemo.outcome,surgery.outcome,observation.outcome)
    value <- cbind(optimal.value,tree.value,guideline2014.value,guideline2018.value,chemo.value,surgery.value,observation.value)
    save(loss,outcome,value,file = filename)
  } else{
    loss[,'cost'] <- as.numeric(loss[,'cost'])
    loss[,'mean'] <- as.numeric(loss[,'mean'])
    loss[,'lower'] <- as.numeric(loss[,'lower'])
    loss[,'upper'] <- as.numeric(loss[,'upper'])
    loss[,'standard_error'] <- as.numeric(loss[,'standard_error'])
    loss[,'QL'] <- as.numeric(loss[,'QL'])
    loss[,'QU'] <- as.numeric(loss[,'QU'])

    outcome[,'cost'] <- as.numeric(outcome[,'cost'])
    outcome[,'mean'] <- as.numeric(outcome[,'mean'])
    outcome[,'lower'] <- as.numeric(outcome[,'lower'])
    outcome[,'upper'] <- as.numeric(outcome[,'upper'])
    outcome[,'standard_error'] <- as.numeric(outcome[,'standard_error'])
    outcome[,'QL'] <- as.numeric(outcome[,'QL'])
    outcome[,'QU'] <- as.numeric(outcome[,'QU'])

    value[,'cost'] <- as.numeric(value[,'cost'])
    value[,'mean'] <- as.numeric(value[,'mean'])
    value[,'lower'] <- as.numeric(value[,'lower'])
    value[,'upper'] <- as.numeric(value[,'upper'])
    value[,'standard_error'] <- as.numeric(value[,'standard_error'])
    value[,'QL'] <- as.numeric(value[,'QL'])
    value[,'QU'] <- as.numeric(value[,'QU'])

    save(loss,outcome,value,file=filename)
  }
}

# return plot of outcome; numeric detail of outcome
# cost_seq for print
PlotPrintLossOutcomeValue <- function(treatmentpercent,cost_seq,filename){
  load(filename)

  loss <- filter(loss, cost<=0.4)
  outcome <- filter(outcome, cost<=0.4)
  value <- filter(value, cost<=0.4)

  loss$estimator <- ordered(loss$estimator, levels = c("chemo","surgery","guideline2018","optimal","tree"))
  outcome$estimator <- ordered(outcome$estimator, levels = c("chemo","surgery","guideline2018","optimal","tree"))
  value$estimator <- ordered(value$estimator, levels = c("chemo","surgery","guideline2018","optimal","tree"))

  p1 <- ggplot(data = loss,aes(x=cost,y=mean,color=estimator))+
    geom_line() +
    geom_smooth(aes(ymin = lower, ymax = upper),stat = "identity") +
    ylab("loss") +
    xlab('tradeoff')

  p2 <- ggplot(data = outcome,aes(x=cost,y=mean,color=estimator))+
    geom_line() +
    geom_smooth(aes(ymin = lower, ymax = upper),stat = "identity") +
    ylab("outcome") +
    xlab('tradeoff')

  #p3 <- ggplot(data = value,aes(x=cost,y=mean,color=estimator))+
  #  geom_line() +
  #  geom_smooth(aes(ymin = lower, ymax = upper),stat = "identity") +
  #  ylab("value")

  p <- ggarrange(p1,p2,nrow = 1,common.legend = TRUE,legend = "right")


  loss <- loss[is.element(loss$cost,cost_seq), ]
  outcome <- outcome[is.element(outcome$cost,cost_seq), ]
  estimator_tpercent_ID_pair <- cbind(match(loss$cost,treatmentpercent$cost),match(loss$estimator,colnames(treatmentpercent)))
  et <- treatmentpercent[estimator_tpercent_ID_pair]

  print.lov <- cbind(as.character(loss$estimator),'&',et,'&',
                     round(loss$mean,digits = 3),'(',round(loss$lower,digits = 3),',',round(loss$upper,digits = 3),')&',
                     round(outcome$mean,digits = 3),'(',round(outcome$lower,digits = 3),',',round(outcome$upper,digits = 3),')\\')
  print.lov <- as.data.frame(print.lov)
  return(list(p,print.lov))
}

PrintTreatmentPercent <- function(dat,tree,cost_seq){
  tp <- NULL
  TP <- NULL
  for (i in cost_seq) {
    TPFT <- round(sum(tree[(tree$trt.effect > i)&(is.na(tree$score)),'size'])/dim(dat)[1],digits = 2)
    TPFO <- round(mean(dat$y>i),digits = 2)
    guideline2018 <- ifelse(dat$guideline2018==factor('unknown'),1,0)
    TPFG <- round(mean(guideline2018),digits = 2)
    tp <- c(i,TPFO,TPFT,TPFG,1,0)
    TP <- rbind(TP,tp)
  }
  colnames(TP) <- c('cost','optimal','tree','guideline2018','chemo','surgery')
  TP <- as.data.frame(TP)
  return(TP)
}

PlotPercentageTreatment <- function(TP,dat){
  PT.guideline <- TP$guideline2018[1]
  #PT.observe <- mean(dat$trt)
  colors <- c("guideline2018" = "purple", "optimal" = "red", "tree" = "blue",'chemo'='orange','surgery'='green')

  p <- ggplot(data=TP,aes(x=cost,y=optimal,color='optimal'))+
    geom_line() +
    geom_line(data=TP,aes(x=cost,y=tree,color='tree')) +
    geom_hline(aes(yintercept = PT.guideline,color='guideline2018'),linetype="dashed")+
    geom_hline(aes(yintercept = 1,color='chemo'),linetype="dashed") +
    geom_hline(aes(yintercept = 0,color='surgery'),linetype="dashed") +
    labs(x = "tradeoff",
         y = "Prop. chemo",
         color = "estimator") +
    scale_color_manual(values = colors)
    return(p)
}



PlotCostSensTreePlain <- function(tree,cost,dat){
  library(partykit)
  library(ggparty)

  tree.pruned <- tree
  tree.pruned$trt <- 1*(tree$trt.effect > cost)
  tree.pruned <- prune_search(1,tree.pruned)

  n <- 0
  ss <- NULL
  # generate ss
  # sss is the partynode for party
  sss <- gsub(',kids=list{}','',ss,fixed = TRUE)
  sss <- gsub('{','(',sss,fixed = TRUE)
  sss <- gsub('}',')',sss,fixed = TRUE)
  # split_party is the partysplit for party
  split_party <- generate_split(dat,tree.pruned)
  eval(parse(text=split_party))
  pn <- eval(parse(text=sss))
  # party node combined
  #py <- party(pn,dat)
  py
  plot(py)

  # ggparty(py,horizontal = FALSE,terminal_space = 0.3,
  #         #layout = data.frame(id=1,x=0.25,y=1),#for cost=0.25
  #         add_vars = list(t =
  #                           function(data, node) {
  #                             list(round(mean(node$data$trt)*100,1))
  #                           },
  #                         ate =
  #                           function(data,node){
  #                             list(round(mean(node$data$y),3))
  #                           }
  #         )
  # ) +
  #   geom_edge() +
  #   geom_edge_label() +
  #   geom_node_label(aes(label = paste0(splitvar,"\n N=", nodesize,'(',t,'%)',' \n ate=',ate)),
  #                   ids = "inner"
  #   ) +
  #   geom_node_label(aes(label = paste0("N=", nodesize,'(',t,'%)',' \n ate=',ate),col=factor(ate<cost)),
  #                   ids = "terminal", nudge_y = -0.33, nudge_x = 0.01
  #   ) +
  #   geom_node_plot(
  #     ids= 'terminal',
  #     gglist = list(geom_density(aes(x=y),adjust=5),
  #                   geom_vline(xintercept = cost,col='red'),
  #                   ylab(''),
  #                   xlab('')
  #     )
  #   )+
  #   theme(legend.position = 'none')

}

SigmaToRisk <- function(R.ave){
  R.ave[3,] <- 1/2*log(2*pi)+R.ave[3,]/2+1/2
  return(R.ave)
}

# c1=(L00^1,L01^1,L10^1,L11^1); c0=(L00^0,L01^0,L10^0,L11^0)
# in marginal loss function parameterization, ct <- 0.05; cd <- 1;
OptimalProb <- function(p1_,p_1,rho,c1,c0){
  p11 <- rho*sqrt(p1_*(1-p1_)*p_1*(1-p_1)) + p1_*p_1
  p10 <- p1_ - p11
  p01 <- p_1 - p11
  p00 <- 1- p11 - p10 - p01

  L1 <- c1[1]*p00+c1[2]*p01+c1[3]*p10+c1[4]*p11
  L0 <- c0[1]*p00+c0[2]*p01+c0[3]*p10+c0[4]*p11
  optimalprob <- apply(L1 < L0,2,mean)
  return(optimalprob)
}


# input: processed data,xvars for growing tree, initial large tree
# output: alpha and error of pruned tree
TreeComplexityVSObjectiveFunctionBernoulli <- function(dat,xvars,tree,min.ndsz,max.depth,pt,loc){
  Alpha_TreeSize <- Get_Alpha_TreeSize(tree)
  Alpha <- Alpha_TreeSize[[1]]
  TreeSize <- Alpha_TreeSize[[2]]
  Risk <- NULL
  split.var <- 1:length(xvars);ctg=6:length(xvars); mtry=length(split.var)
  yhat <- mean(dat$optimal.prob)
  SLL <- SumLLBernoulii(dat$optimal.prob,yhat)
  N <- dim(dat)[1]

  for (j in 1:length(Alpha)){
    tree <- grow.INT.ILL.OptimalDecision(dat,
                         split.var,ctg,
                         min.ndsz=min.ndsz, pt=pt, max.depth,
                         mtry,
                         alpha = Alpha[j],loc)
    risk <- (sum(tree$score,na.rm = TRUE) + SLL)/N
    Risk <- c(Risk,risk)
  }
  Risk <- rbind(Alpha,TreeSize,Risk)
  return(Risk)
}

SumLLBernoulii <- function(y,yhat){
  SLL <- sum(y*log(yhat)+(1-y)*log(1-yhat))
  return(SLL)
}

Permute_VarIP_Bernouli <- function(bootloop,dat,xvars,min.ndsz,pt,max.depth){
  ctg <- 6:length(xvars)
  n_sam <- nrow(dat)
  Importance.permute <- NULL
  pb <- txtProgressBar(min = 0, max = bootloop, style = 3)

  set.seed(123)
  for (i in seq(bootloop)) {
    train.id <- sample(1:nrow(dat),n_sam,replace = TRUE)
    train <- dat[train.id,]
    tree <- TreeGrownAndPruneOptimalDecision(train,xvars,min.ndsz=min.ndsz,pt=pt,max.depth=max.depth,alpha=0,loc=FALSE)
    #2. compute J(T_b)
    score.bench <- test.Jscore(train, tree, ctg)
    #3. permute x_{j} then compute J(T_b)_j, and compute J(T_b)-J(T_b)_j
    importance.permute <- NULL
    for (xvar in xvars) {
      permute.var <- permute(train[,xvar])
      test.permute <- train
      test.permute[,xvar] <- permute.var
      score.permute <- test.Jscore(test.permute,tree,ctg)
      importance.var <- (score.bench - score.permute)/abs(score.bench)
      importance.permute <- c(importance.permute,importance.var)
    }
    Importance.permute <- rbind(Importance.permute,importance.permute)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(Importance.permute)
}

test.Jscore <- function(test, tree,ctg){
  N <- dim(test)[1]
  Outcome <- predict_initree_OptimalDecision(tree,test,ctg)
  Risk.test <- sum(test$optimal.prob*log(Outcome)+(1-test$optimal.prob)*log(1-Outcome))/N
  return(Risk.test)
}

PlotFI_Bernoulli <- function(bootstrap.IP,xvars){
  bootstrap.IP <- as.data.frame(bootstrap.IP)
  bootstrap.IP <- apply(bootstrap.IP,2, mean)
  bootstrap.IP <- as.data.frame(bootstrap.IP)
  bootstrap.IP<-  as.data.frame(cbind(xvars,bootstrap.IP))
  colnames(bootstrap.IP) <- c("var",'importance')
  # convert factor to numeric: factor->character->numeric
  bootstrap.IP$importance <- as.numeric(as.character(bootstrap.IP$importance))
  bootstrap.IP <- filter(bootstrap.IP,importance!=0)
  # order var based on importance ranking
  bootstrap.IP$var <- with(bootstrap.IP, factor(var,
                                                levels=var[order(importance)]))
  pFI <- ggplot(data=bootstrap.IP, aes(x=var, y=importance)) +
    geom_bar(stat="identity") + coord_flip() +
    theme(axis.title.y=element_blank(),
          text = element_text(size=10),
          axis.title.x = element_text(size = 18))
  #ggtitle("Feature importance via random forest")

  return(pFI)
}

Compute_loss_value_Bernoulli <- function(year,data,xvars,min.ndsz,pt,max.depth,alpha,p1_,p_1,rho,c1,c0,ct_seq){
  filename <- paste('loss_outcome_rho',rho,'_',year,'fp.RData',sep='')
  ctg <- 6:length(xvars)
  # =============================================================
  # correlation of Y(1),Y(0): rho
  # theta
  # =============================================================
  rho <- rho
  p11 <- rho*sqrt(p1_*(1-p1_)*p_1*(1-p_1)) + p1_*p_1
  p10 <- p1_ - p11
  p01 <- p_1 - p11
  p00 <- 1- p11 - p10 - p01

  expected.y1 <- apply(p1_, 2, mean); var.y1 <- apply(p1_,2,var)
  expected.y0 <- apply(p_1, 2, mean); var.y0 <- apply(p_1,2,var)
  L0 <- c0[1]*p00+c0[2]*p01+c0[3]*p10+c0[4]*p11
  expected.L0 <- apply(L0, 2, mean); var.L0 <- apply(L0, 2, var)

  # =============================================================
  # COMPUTE POSTERIOR EXPECTED LOSS AND POSTERIOR EXPECTED OUTCOME
  # =============================================================
  guideline2018.action <- ifelse(data$guideline2018=='unknown',ifelse(data$y>0,1,0),data$guideline2018)

  # =============================================================
  loss <- data.frame(cost=numeric(),estimator=character(),mean=numeric(),
                     lower=numeric(),upper=numeric(),standard_error=numeric(),
                     QL=numeric(),QU=numeric(),stringsAsFactors = FALSE)
  outcome <- data.frame(cost=numeric(),estimator=character(),mean=numeric(),
                        lower=numeric(),upper=numeric(),standard_error=numeric(),
                        QL=numeric(),QU=numeric(),stringsAsFactors = FALSE)

  i <- 0
  for (delta in ct_seq){ # from 0 to 1 (include 0 an 1)
    c1_update <- c1+delta
    optimalprob <- OptimalProb(p1_,p_1,rho,c1_update,c0)
    data$optimal.prob <- optimalprob
    tree <- TreeGrownAndPruneOptimalDecision(data,xvars,min.ndsz=min.ndsz,pt=pt,max.depth=max.depth,alpha=alpha,loc=FALSE)
    tree.decision <- predict_initree_OptimalDecision(tree,data,ctg)

    threshold <- 0.5
    L1 <- c1_update[1]*p00+c1_update[2]*p01+c1_update[3]*p10+c1_update[4]*p11;
    expected.L1 <- apply(L1, 2, mean); var.L1 <- apply(L1, 2, var)

    chemo.outcome <- expected.y1
    chemo.loss <- expected.L1
    chemo.outcome.var <- var.y1
    chemo.loss.var <- var.L1

    surgery.outcome <- expected.y0
    surgery.loss <- expected.L0
    surgery.outcome.var <- var.y0
    surgery.loss.var <- var.L0

    guideline2018.outcome <- ifelse(guideline2018.action==1,expected.y1,expected.y0)
    guideline2018.loss <- ifelse(guideline2018.action==1,expected.L1,expected.L0)
    guideline2018.outcome.var <- ifelse(guideline2018.action==1,var.y1,var.y0)
    guideline2018.loss.var <- ifelse(guideline2018.action==1,var.L1,var.L0)

    optimal.action <- ifelse(data$optimal.prob>threshold,1,0)
    optimal.outcome <- ifelse(optimal.action==1,expected.y1,expected.y0)
    optimal.loss <- ifelse(optimal.action==1,expected.L1,expected.L0)
    optimal.outcome.var <- ifelse(optimal.action==1,var.y1,var.y0)
    optimal.loss.var <- ifelse(optimal.action==1,var.L1,var.L0)

    tree.action <- ifelse(tree.decision>threshold,1,0)
    tree.outcome <- ifelse(tree.action==1,expected.y1,expected.y0)
    tree.loss <- ifelse(tree.action==1,expected.L1,expected.L0)
    tree.outcome.var <- ifelse(tree.action==1,var.y1,var.y0)
    tree.loss.var <- ifelse(tree.action==1,var.L1,var.L0)


    i <- i+1
    loss[i,] <- append.data(chemo.loss,'chemo',delta,chemo.loss.var)
    outcome[i,] <- append.data(chemo.outcome,'chemo',delta,chemo.outcome.var)

    i <- i+1
    loss[i,] <- append.data(surgery.loss,'surgery',delta,surgery.loss.var)
    outcome[i,] <- append.data(surgery.outcome,'surgery',delta,surgery.outcome.var)

    i <- i+1
    loss[i,] <- append.data(guideline2018.loss,'guideline2018',delta,guideline2018.loss.var)
    outcome[i,] <- append.data(guideline2018.outcome,'guideline2018',delta,guideline2018.outcome.var)

    i <- i+1
    loss[i,] <- append.data(optimal.loss,'optimal',delta,optimal.loss.var)
    outcome[i,] <- append.data(optimal.outcome,'optimal',delta,optimal.outcome.var)

    i <- i+1
    loss[i,] <- append.data(tree.loss,'tree',delta,tree.loss.var)
    outcome[i,] <- append.data(tree.outcome,'tree',delta,tree.outcome.var)

    print(delta)
  }

  if(length(ct_seq)==1){
    loss <- cbind(optimal.loss,tree.loss,guideline2018.loss,chemo.loss,surgery.loss,observation.loss)
    outcome <- cbind(optimal.outcome,tree.outcome,guideline2018.outcome,chemo.outcome,surgery.outcome,observation.outcome)
    save(loss,outcome,file = filename)
  } else{
    loss[,'cost'] <- as.numeric(loss[,'cost'])
    loss[,'mean'] <- as.numeric(loss[,'mean'])
    loss[,'lower'] <- as.numeric(loss[,'lower'])
    loss[,'upper'] <- as.numeric(loss[,'upper'])
    loss[,'standard_error'] <- as.numeric(loss[,'standard_error'])
    loss[,'QL'] <- as.numeric(loss[,'QL'])
    loss[,'QU'] <- as.numeric(loss[,'QU'])

    outcome[,'cost'] <- as.numeric(outcome[,'cost'])
    outcome[,'mean'] <- as.numeric(outcome[,'mean'])
    outcome[,'lower'] <- as.numeric(outcome[,'lower'])
    outcome[,'upper'] <- as.numeric(outcome[,'upper'])
    outcome[,'standard_error'] <- as.numeric(outcome[,'standard_error'])
    outcome[,'QL'] <- as.numeric(outcome[,'QL'])
    outcome[,'QU'] <- as.numeric(outcome[,'QU'])

    save(loss,outcome,file=filename)
  }
}

ComputePercentTreatment_Bernoulli <- function(p1_,p_1,rho,ct_seq,alpha,data,xvars,min.ndsz,pt,max.depth){
  tp <- NULL
  TP <- NULL
  ctg <- 6:length(xvars)

  for (delta in ct_seq) {
    c1_update <- c1+delta
    optimalprob <- OptimalProb(p1_,p_1,rho,c1_update,c0)
    data$optimal.prob <- optimalprob
    tree <- TreeGrownAndPruneOptimalDecision(data,xvars,min.ndsz=min.ndsz,pt=pt,max.depth=max.depth,alpha=alpha,loc=FALSE)
    tree.decision <- predict_initree_OptimalDecision(tree,data,ctg)

    TPFT <- round(mean(tree.decision>0.5),digits = 2)
    TPFO <- round(mean(data$optimal.prob>0.5),digits = 2)
    guideline2018 <- ifelse(dat$guideline2018==factor('unknown'),1,0)
    TPFG <- round(mean(guideline2018),digits = 2)
    tp <- c(delta,TPFO,TPFT,TPFG,1,0)
    TP <- rbind(TP,tp)
  }
  colnames(TP) <- c('cost','optimal','tree','guideline2018','chemo','surgery')
  TP <- as.data.frame(TP)
  return(TP)
}

PlotLossOutcome_Bernoulli <- function(treatmentpercent,cost_seq,filename){
  load(filename)

  loss$estimator <- ordered(loss$estimator, levels = c("chemo","surgery","guideline2018","optimal","tree"))
  outcome$estimator <- ordered(outcome$estimator, levels = c("chemo","surgery","guideline2018","optimal","tree"))

  p1 <- ggplot(data = loss,aes(x=cost,y=mean,color=estimator))+
    geom_line() +
    geom_smooth(aes(ymin = lower, ymax = upper),stat = "identity") +
    ylab("loss") +
    xlab('tradeoff') +
    theme(legend.title=element_blank(),
          text = element_text(size=18))

  p2 <- ggplot(data = outcome,aes(x=cost,y=mean,color=estimator))+
    geom_line() +
    geom_smooth(aes(ymin = lower, ymax = upper),stat = "identity") +
    ylab("value") +
    xlab('tradeoff') +
    theme(legend.title=element_blank(),
          text = element_text(size=18))

  #p <- ggarrange(p1,p2,nrow = 1,common.legend = TRUE,legend = "right")

  percent <- as.data.frame(NULL)
  for(i in cost_seq){
    cost_tpercent <- t(treatmentpercent[treatmentpercent$cost==i,c('chemo','surgery','guideline2018','optimal','tree')])
    percent <- c(percent,cost_tpercent)
  }

  print.lov <- cbind(as.character(loss$estimator),'&',percent,'&',
                     round(loss$mean,digits = 3),'(',round(loss$standard_error*100,digits = 3),')&',
                     round(outcome$mean,digits = 3),'(',round(outcome$standard_error*100,digits = 3),')\\')
  print.lov <- as.data.frame(print.lov)
  colnames(print.lov) <- c('decision rules','','percent','','loss mean','','var','','outcome mean','','var','')
  rownames(print.lov) <- NULL
  print(print.lov,row.names = FALSE)

  return(list(p1,p2))
}

Plot4PaperTreePolicy <- function(tree,dat){
  library(partykit)
  library(ggparty)
  tree.pruning <- tree
  n <- 0
  ss <- NULL
  tree_search(1,tree.pruning) # generate ss
  # sss is the partynode for party
  sss <- gsub(',kids=list{}','',ss,fixed = TRUE)
  sss <- gsub('{','(',sss,fixed = TRUE)
  sss <- gsub('}',')',sss,fixed = TRUE)
  # split_party is the partysplit for party
  split_party <- generate_split(dat,tree.pruning)
  eval(parse(text=split_party))
  pn <- eval(parse(text=sss))
  # party node combined
  py <- party(pn,dat)



  P.tree <- ggparty(py,horizontal = FALSE,terminal_space = 0.3,
              #layout = data.frame(id=1,x=0.25,y=1),#for cost=0.25
              add_vars = list(t =
                                function(data, node) {
                                  list(round(mean(node$data$trt)*100,1))
                                },
                              p.tree =
                                function(data,node){
                                  list(round(mean(node$data$optimal.prob),3))
                                }
              )
      ) +
        geom_edge() +
        geom_edge_label() +
        geom_node_label(aes(label = paste(splitvar,"\n N=", nodesize, "\n",p.tree,sep="")),
                        ids = "inner"
        ) +
        geom_node_label(aes(label = paste("N=", nodesize, "\n",p.tree),col=factor(p.tree<0.5)),
                        ids = "terminal", nudge_y = -0.33, nudge_x = 0.01
        ) +
        geom_node_plot(
          ids= 'terminal',
          gglist = list(geom_density(aes(x=optimal.prob),adjust=5),
                        geom_vline(xintercept = 0.5,col='red'),
                        ylim(0,5),
                        ylab(''),
                        xlab('')
          )
        )+
        theme(legend.position = 'none')
  P.tree

  return(p.tree)
}
#paste(splitvar,"\n N=", nodesize, "\n p(z^t=1)=",p.tree,sep=""))

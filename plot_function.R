d.m23.acc <- function(complexity, d.m.acc, d.m.estime, d.m.true){
  n <- dim(d.m.true)[1]
  for(i in seq(complexity)){
    d.m.estime.i <- t(d.m.estime[i,,])
    d.m.mean.4sim <- apply((d.m.estime.i == d.m.true)*1, 1, mean)
    mean <- mean(d.m.mean.4sim)
    se <- sd(d.m.mean.4sim)
    d.m.acc[i,"complexity"] <- i
    d.m.acc[i,"mean"] <- mean
    d.m.acc[i,"upper"] <- mean + 1.98*se/sqrt(n)
    d.m.acc[i,"lower"] <- mean - 1.98*se/sqrt(n)
  }
  return(d.m.acc)
}

d.m23.precision <- function(complexity, d.m.acc, d.m.estime, d.m.true){
  for(i in seq(complexity)){
    d.m.estime.i <- t(d.m.estime[i,,])
    d.m.estime.predict <- apply((d.m.estime.i == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.m.true) & (d.m.true==1))*1,1,sum)
    precision <- d.m.estime.true/d.m.estime.predict
    mean <- mean(na.omit(precision))
    sd <- sd(na.omit(precision))
    d.m.acc[i,"complexity"] <- i
    d.m.acc[i,"mean"] <- mean
    d.m.acc[i,"upper"] <- mean + 1.98*sd/sqrt(sum(!is.na(precision)))
    d.m.acc[i,"lower"] <- mean - 1.98*sd/sqrt(sum(!is.na(precision)))
  }
  return(d.m.acc)
}

d.m23.recall <- function(complexity, d.m.acc, d.m.estime, d.m.true){
  for(i in seq(complexity)){
    d.m.estime.i <- t(d.m.estime[i,,])
    d.m.true.positive <- apply((d.m.true == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.m.true) & (d.m.true==1))*1,1,sum)
    recall <- d.m.estime.true/d.m.true.positive
    mean <- mean(na.omit(recall))
    sd <- sd(na.omit(recall))
    d.m.acc[i,"complexity"] <- i
    d.m.acc[i,"mean"] <- mean
    d.m.acc[i,"upper"] <- mean + 1.98*sd/sum(!is.na(recall))
    d.m.acc[i,"lower"] <- mean - 1.98*sd/sum(!is.na(recall))
  }
  return(d.m.acc)
}

d.m23.f1 <- function(complexity, d.m.acc, d.m.estime, d.m.true){
  #browser()
  for(i in seq(complexity)){
    d.m.estime.i <- t(d.m.estime[i,,])
    d.m.estime.predict <- apply((d.m.estime.i == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.m.true) & (d.m.true==1))*1,1,sum)
    d.m.precision <- d.m.estime.true/d.m.estime.predict
    
    d.m.true.positive <- apply((d.m.true == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.m.true) & (d.m.true==1))*1,1,sum)
    d.m.recall <- d.m.estime.true/d.m.true.positive
    
    d.m.f1 <- 2*d.m.precision*d.m.recall/(d.m.precision+d.m.recall)
    mean <- mean(na.omit(d.m.f1))
    sd <- sd(na.omit(d.m.f1))
    
    d.m.acc[i,"complexity"] <- i
    d.m.acc[i,"mean"] <- mean
    d.m.acc[i,"upper"] <- mean + 1.98*sd/sqrt(sum(!is.na(d.m.f1)))
    d.m.acc[i,"lower"] <- mean - 1.98*sd/sqrt(sum(!is.na(d.m.f1)))
  }
  return(d.m.acc)
}

d.loss <- function(complexity,d.m.estime,m.loss,L1,L0){
  n <- dim(d.m.estime)[3]
  for (i in seq(complexity)) {
    d.m.estime.i <- t(d.m.estime[i,,])
    m.L <- apply(L1*d.m.estime.i + L0*(1-d.m.estime.i),1,mean)
    mean <- mean(m.L)
    sd <- sd(m.L)
    m.loss[i,"complexity"] <- i
    m.loss[i,"mean"] <- mean
    m.loss[i,"upper"] <- mean + 1.98*sd/sqrt(n)
    m.loss[i,"lower"] <- mean - 1.98*sd/sqrt(n)
  }
  return(m.loss)
}

d.m23outcome <- function(complexity, m.outcome, d.m.estime, y1, y0){
  #browser()
  n <- dim(d.m.estime)[3]
  for (i in seq(complexity)) {
    d.m.estime.i <- t(d.m.estime[i,,])
    m.outcome.i <- apply(y1*d.m.estime.i + y0*(1-d.m.estime.i),1,mean)
    mean <- mean(m.outcome.i)
    sd <- sd(m.outcome.i)
    m.outcome[i,"complexity"] <- i
    m.outcome[i,"mean"] <- mean
    m.outcome[i,"upper"] <- mean + 1.98*sd/sqrt(n)
    m.outcome[i,"lower"] <- mean - 1.98*sd/sqrt(n)
  }
  return(m.outcome)
}


plot_loss_outcome <- function(outcome,complexity,threhold){
  # d.true is an array with 2 dimension, n_simulation * n_sample,
  # d.m2 and d.m3 are arrays with 3 dimension, complexity*n_sample*n_simulation
  d.true <- (outcome[[1]]-outcome[[2]] > threhold)*1
  d.m1 <- (outcome[[3]]-outcome[[4]] > threhold)*1 
  d.m2 <- (outcome[[5]] > threhold)*1
  d.m3 <- (outcome[[6]]-outcome[[7]] > threhold)*1
  
  m2.loss <- m3.loss <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(m2.loss) <- c("complexity","mean","upper","lower")
  colnames(m3.loss) <- c("complexity","mean","upper","lower")
  p11 <- 0*sqrt(outcome[[1]]*(1-outcome[[1]])*outcome[[2]]*(1-outcome[[2]])) + outcome[[1]]*outcome[[2]]
  p10 <- outcome[[1]] - p11
  p01 <- outcome[[2]] - p11
  p00 <- 1 - p11 - p10 - p01
  L1 <- p00*(1+threhold) + p01*(1+threhold) +p10*threhold + p11*threhold
  L0 <- p00 + p10
  
  m1.loss <- apply(L1*d.m1 + L0*(1-d.m1),1,mean)
  m1.loss.mean <- mean(m1.loss)
  m1.loss.sd <- sd(m1.loss)
  m1.loss.upper <- m1.loss.mean + 1.98*m1.loss.sd/sqrt(100)
  m1.loss.lower <- m1.loss.mean - 1.98*m1.loss.sd/sqrt(100)
  m2.loss <- d.loss(complexity,d.m2, m2.loss, L1, L0)
  m3.loss <- d.loss(complexity,d.m3, m3.loss, L1, L0)
  
  m2.outcome <- m3.outcome <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(m2.outcome) <- c("complexity","mean","upper","lower")
  colnames(m3.outcome) <- c("complexity","mean","upper","lower")
  m1.outcome <- apply(outcome[[1]]*d.m1 + outcome[[2]]*(1-d.m1),1,mean)
  m1.outcome.mean <- mean(m1.outcome)
  m1.outcome.sd <- sd(m1.outcome)
  m1.outcome.upper <- m1.outcome.mean + 1.98*m1.outcome.sd/sqrt(100)
  m1.outcome.lower <- m1.outcome.mean + 1.98*m1.outcome.sd/sqrt(100)
  m2.outcome <- d.m23outcome(complexity, m2.outcome, d.m2, y1=outcome[[1]], y0=outcome[[2]])
  m3.outcome <- d.m23outcome(complexity, m3.outcome, d.m3, y1=outcome[[1]], y0=outcome[[2]])
  
  colors <- c("f1" = "black","f2" = "blue", "f3" = "red")
  shapes <- c("f1" = 1, "f2" = 0, "f3" = 2)
  
  p.loss <- ggplot() +
    geom_ribbon(aes(x=seq(complexity),y=m1.loss.mean,ymin=m1.loss.lower, ymax=m1.loss.upper),
                fill = "grey70",alpha=0.8)+
    geom_point(aes(x=seq(complexity),y=m2.loss$mean,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=m3.loss$mean,color="f3"),size=3)+
    geom_line(aes(x=seq(complexity),y=m2.loss$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=m3.loss$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_hline(aes(yintercept=m1.loss.mean, color="f1"),linetype = "dashed") +
    geom_errorbar(aes(x=seq(complexity), ymin=m2.loss$lower, ymax=m2.loss$upper, color="f2"), width=0.5, size=0.3)+
    geom_errorbar(aes(x=seq(complexity), ymin=m3.loss$lower, ymax=m3.loss$upper, color="f3"), width=0.5, size=0.3)+
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{"sim"}))+
    labs(x = 'depth',
         y = 'loss',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))
  
  p.outcome <- ggplot() +
    geom_ribbon(aes(x=seq(complexity),y=m1.outcome.mean,ymin=m1.outcome.lower, ymax=m1.outcome.upper),
                fill = "grey70",alpha=0.8)+
    geom_point(aes(x=seq(complexity),y=m2.outcome$mean,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=m3.outcome$mean,color="f3"),size=3)+
    geom_line(aes(x=seq(complexity),y=m2.outcome$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=m3.outcome$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_hline(aes(yintercept=m1.outcome.mean, color="f1"),linetype = "dashed") +
    geom_errorbar(aes(x=seq(complexity), ymin=m2.outcome$lower, ymax=m2.outcome$upper, color="f2"), width=0.5, size=0.3)+
    geom_errorbar(aes(x=seq(complexity), ymin=m3.outcome$lower, ymax=m3.outcome$upper, color="f3"), width=0.5, size=0.3)+
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{"sim"}))+
    labs(x = 'depth',
         y = 'outcome',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))
  
  return(list(p.loss,p.outcome))
}

# return accuracy, F1, precision, recall
plot_accuracy_metrics <- function(outcome,complexity,threhold){
  #browser()
  d.true <- (outcome[[1]]-outcome[[2]] > threhold)*1
  d.m1 <- (outcome[[3]]-outcome[[4]] > threhold)*1 
  d.m2 <- (outcome[[5]] > threhold)*1
  d.m3 <- (outcome[[6]]-outcome[[7]] > threhold)*1
  
  d.m1.acc.mean <- mean(apply((d.m1==d.true)*1, 1, mean))
  d.m1.acc.sd <- sd(apply((d.m1==d.true)*1, 1, mean))
  d.m1.acc.upper <- d.m1.acc.mean + 1.98*d.m1.acc.sd/sqrt(100)
  d.m1.acc.lower <- d.m1.acc.mean - 1.98*d.m1.acc.sd/sqrt(100)

  d.m2.acc <- d.m3.acc <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(d.m2.acc) <- c("complexity","mean","upper","lower")
  colnames(d.m3.acc) <- c("complexity","mean","upper","lower")
  d.m2.acc <- d.m23.acc(complexity = 10, d.m.acc = d.m2.acc, d.m.estime = d.m2, d.m.true = d.true)
  d.m3.acc <- d.m23.acc(complexity = 10, d.m.acc = d.m3.acc, d.m.estime = d.m3, d.m.true = d.true)

  d.m1.true.positive <- apply(((d.m1==d.true)&(d.true==1))*1,1,sum)
  d.m1.positive <- apply(d.m1,1,sum)
  d.m1.precision <- d.m1.true.positive/d.m1.positive
  d.m1.precision.mean <- mean(d.m1.precision)
  d.m1.precision.sd <- sd(d.m1.precision)
  d.m1.precision.upper <- d.m1.precision.mean+1.98*d.m1.precision.sd/sqrt(100)
  d.m1.precision.lower <- d.m1.precision.mean-1.98*d.m1.precision.sd/sqrt(100)

  d.m2.precision <- d.m3.precision <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(d.m2.precision) <- c("complexity","mean","upper","lower")
  colnames(d.m3.precision) <- c("complexity","mean","upper","lower")
  d.m2.precision <- d.m23.precision(complexity = 10, d.m.acc = d.m2.precision, d.m.estime = d.m2, d.m.true = d.true)
  d.m3.precision <- d.m23.precision(complexity = 10, d.m.acc = d.m3.precision, d.m.estime = d.m3, d.m.true = d.true)

  d.m1.true.positive <- apply(((d.m1==d.true)&(d.true==1))*1,1,sum)
  d.true.positive <- apply(d.true,1,sum)
  d.m1.recall <- d.m1.true.positive/d.true.positive
  d.m1.recall.mean <- mean(d.m1.recall)
  d.m1.recall.sd <- sd(d.m1.recall)
  d.m1.recall.upper <- d.m1.recall.mean+1.98*d.m1.recall.sd/sqrt(100)
  d.m1.recall.lower <- d.m1.recall.mean-1.98*d.m1.recall.sd/sqrt(100)
  #
  d.m2.recall <- d.m3.recall <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(d.m2.recall) <- c("complexity","mean","upper","lower")
  colnames(d.m3.recall) <- c("complexity","mean","upper","lower")
  d.m2.recall <- d.m23.recall(complexity = 10, d.m.acc = d.m2.recall, d.m.estime = d.m2, d.m.true = d.true)
  d.m3.recall <- d.m23.recall(complexity = 10, d.m.acc = d.m3.recall, d.m.estime = d.m3, d.m.true = d.true)


  d.m1.f1 <- 2*d.m1.precision*d.m1.recall/(d.m1.precision+d.m1.recall)
  d.m1.f1.mean <- mean(d.m1.f1)
  d.m1.f1.sd <- sd(d.m1.f1)
  d.m1.f1.upper <- d.m1.f1.mean+1.98*d.m1.f1.sd/sqrt(100)
  d.m1.f1.lower <- d.m1.f1.mean-1.98*d.m1.f1.sd/sqrt(100)
  
  d.m2.f1 <- d.m3.f1 <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(d.m2.recall) <- c("complexity","mean","upper","lower")
  colnames(d.m3.recall) <- c("complexity","mean","upper","lower")
  d.m2.f1 <- d.m23.f1(complexity = 10, d.m.acc = d.m2.f1, d.m.estime = d.m2, d.m.true = d.true)
  d.m3.f1 <- d.m23.f1(complexity = 10, d.m.acc = d.m3.f1, d.m.estime = d.m3, d.m.true = d.true)
  
  colors <- c("f1" = "black","f2" = "blue", "f3" = "red")
  shapes <- c("f1" = 1, "f2" = 0, "f3" = 2)
  
  p.acc <- ggplot() +
    geom_ribbon(aes(x=seq(complexity),y=d.m1.acc.mean,ymin=d.m1.acc.lower, ymax=d.m1.acc.upper),
                fill = "grey70",alpha=0.8)+
    geom_point(aes(x=seq(complexity),y=d.m2.acc$mean,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=d.m3.acc$mean,color="f3"),size=3)+
    geom_line(aes(x=seq(complexity),y=d.m2.acc$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=d.m3.acc$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_hline(aes(yintercept=d.m1.acc.mean, color="f1"),linetype = "dashed") +
    geom_errorbar(aes(x=seq(complexity), ymin=d.m2.acc$lower, ymax=d.m2.acc$upper, color="f2"), width=1, size=0.5)+
    geom_errorbar(aes(x=seq(complexity), ymin=d.m3.acc$lower, ymax=d.m3.acc$upper, color="f3"), width=1, size=0.5)+
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{'sim'}))+
    labs(x = 'depth',
         y = 'acc.mean',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))
  
  p.precision <- ggplot() +
    geom_ribbon(aes(x=seq(complexity),y=d.m1.precision.mean,ymin=d.m1.precision.lower, ymax=d.m1.precision.upper),
                fill = "grey70",alpha=0.8)+
    geom_point(aes(x=seq(complexity),y=d.m2.precision$mean,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=d.m3.precision$mean,color="f3"),size=3)+
    geom_line(aes(x=seq(complexity),y=d.m2.precision$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=d.m3.precision$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_hline(aes(yintercept=d.m1.precision.mean, color="f1"),linetype = "dashed") +
    geom_errorbar(aes(x=seq(complexity), ymin=d.m2.precision$lower, ymax=d.m2.precision$upper, color="f2"), width=1, size=0.5)+
    geom_errorbar(aes(x=seq(complexity), ymin=d.m3.precision$lower, ymax=d.m3.precision$upper, color="f3"), width=1, size=0.5)+
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{'sim'}))+
    labs(x = 'depth',
         y = 'precision.mean',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))

  p.recall <- ggplot() +
    geom_ribbon(aes(x=seq(complexity),y=d.m1.recall.mean,ymin=d.m1.recall.lower, ymax=d.m1.recall.upper),
                fill = "grey70",alpha=0.8)+
    geom_point(aes(x=seq(complexity),y=d.m2.recall$mean,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=d.m3.recall$mean,color="f3"),size=3)+
    geom_line(aes(x=seq(complexity),y=d.m2.recall$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=d.m3.recall$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_hline(aes(yintercept=d.m1.recall.mean, color="f1"),linetype = "dashed") +
    geom_errorbar(aes(x=seq(complexity), ymin=d.m2.recall$lower, ymax=d.m2.recall$upper, color="f2"), width=1, size=0.5)+
    geom_errorbar(aes(x=seq(complexity), ymin=d.m3.recall$lower, ymax=d.m3.recall$upper, color="f3"), width=1, size=0.5)+
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{'sim'}))+
    labs(x = 'depth',
         y = 'recall.mean',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))


  p.f1 <- ggplot() +
    geom_ribbon(aes(x=seq(complexity),y=d.m1.f1.mean,ymin=d.m1.f1.lower, ymax=d.m1.f1.upper),
                fill = "grey70",alpha=0.8)+
    geom_point(aes(x=seq(complexity),y=d.m2.f1$mean,color="f2"),size=3)+
    geom_point(aes(x=seq(complexity),y=d.m3.f1$mean,color="f3"),size=3)+
    geom_line(aes(x=seq(complexity),y=d.m2.f1$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=seq(complexity),y=d.m3.f1$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_hline(aes(yintercept=d.m1.f1.mean, color="f1"),linetype = "dashed") +
    geom_errorbar(aes(x=seq(complexity), ymin=d.m2.f1$lower, ymax=d.m2.f1$upper, color="f2"), width=1, size=0.5)+
    geom_errorbar(aes(x=seq(complexity), ymin=d.m3.f1$lower, ymax=d.m3.f1$upper, color="f3"), width=1, size=0.5)+
    scale_x_continuous(breaks = seq(complexity))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{'sim'}))+
    labs(x = 'depth',
         y = 'f1.mean',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))
  
  return(list(p.acc,p.f1,p.precision,p.recall))
}


# =======================================================
# BART-logistic_regression functions
# =======================================================
compute_bivariate_prob_decision <- function(p1_,p_1, 
                                            ct00,ct01,ct10,ct11,
                                            cc00,cc01,cc10,cc11,rho){
  p11 <- rho*sqrt(p1_*(1-p1_)*p_1*(1-p_1)) + p1_*p_1
  p10 <- p1_ - p11
  p01 <- p_1 - p11
  p00 <- 1-p11-p10-p01
  L.1 <- p11*ct11 + p10*ct10+ p01*ct01 + p00*ct00
  L.0 <- p11*cc11 + p10*cc10+ p01*cc01 + p00*cc00
  d.true <- (L.1 < L.0)*1
  return(d.true)
}


d.precision <- function(d.true, d.predict){
  d.predict.1 <- apply((d.predict == 1)*1, 1,sum)
  d.true.1 <- apply(((d.true == d.predict) & (d.true==1))*1,1,sum)
  precision <- d.true.1/d.predict.1
  precision.mean <- mean(na.omit(precision))
  precision.sd <- sd(na.omit(precision))
  precision.upper <- precision.mean + 1.98*precision.sd/sqrt(sum(!is.na(precision)))
  precision.lower <- precision.mean - 1.98*precision.sd/sqrt(sum(!is.na(precision)))
  precision.sd <- round(precision.sd, digits = 3)
  return(c(precision.mean,precision.upper,precision.lower,precision.sd))
}

d.recall <- function(d.true, d.predict){
  d.true.1 <- apply((d.true == 1)*1,1,sum)
  d.predict.1 <- apply(((d.predict == d.true) & (d.true==1))*1,1,sum)
  recall <- d.predict.1/d.true.1
  recall.mean <- mean(na.omit(recall))
  recall.sd <- sd(na.omit(recall))
  recall.upper <- recall.mean + 1.98*recall.sd/sum(!is.na(recall))
  recall.lower <- recall.mean - 1.98*recall.sd/sum(!is.na(recall))
  recall.sd <- round(recall.sd, digits = 3)
  return(c(recall.mean,recall.upper, recall.lower,recall.sd))
}

d.f1 <- function(d.true, d.predict){
  d.predict.1 <- apply((d.predict == 1)*1, 1,sum)
  d.true.1 <- apply(((d.true == d.predict) & (d.true==1))*1,1,sum)
  precision <- d.true.1/d.predict.1
  
  d.true.1 <- apply((d.true == 1)*1,1,sum)
  d.predict.1 <- apply(((d.predict == d.true) & (d.true==1))*1,1,sum)
  recall <- d.predict.1/d.true.1
  
  f1 <- 2*precision*recall/(precision+recall)
  f1.mean <- mean(na.omit(f1))
  f1.sd <- sd(na.omit(f1))
  f1.upper <- f1.mean + 1.98*f1.sd/sqrt(sum(!is.na(f1)))
  f1.lower <- f1.mean - 1.98*f1.sd/sqrt(sum(!is.na(f1)))
  f1.sd <- round(f1.sd, digits = 3)
  return(c(f1.mean,f1.upper,f1.lower, f1.sd))
}

d.loss.2 <- function(d,
                     p1_,p_1,rho,
                     ct00,ct01,ct10,ct11,
                     cc00,cc01,cc10,cc11){
  n <- dim(d)[1]
  p11 <- rho*sqrt(p1_*(1-p1_)*p_1*(1-p_1)) + p1_*p_1
  p10 <- p1_ - p11
  p01 <- p_1 - p11
  p00 <- 1-p11-p10-p01
  L1 <- p11*ct11 + p10*ct10+ p01*ct01 + p00*ct00
  L0 <- p11*cc11 + p10*cc10+ p01*cc01 + p00*cc00
  
  L <- apply(L1*d + L0*(1-d),1,mean)
  L.mean <- round(mean(L),digits = 3)
  L.sd <- sd(L)
  L.lower <- round(L.mean - 1.98*L.sd/sqrt(n),digits = 3)
  L.upper <- round(L.mean + 1.98*L.sd/sqrt(n),digits = 3)
  L.sd <- round(L.sd, digits = 3)
  return(c(L.mean, L.lower, L.upper,L.sd))
}

d.outcome.2 <- function(d,
                        p1_,p_1){
  n <- dim(d)[1]
  O <- apply(p1_*d + p_1*(1-d),1,mean)
  O.mean <- round(mean(O),digits = 3)
  O.sd <- sd(O)
  O.lower <- round(O.mean - 1.98*O.sd/sqrt(n),digits = 3)
  O.upper <- round(O.mean + 1.98*O.sd/sqrt(n),digits = 3)
  O.sd <- round(O.sd, digits = 3)
  return(c(O.mean, O.lower, O.upper, O.sd))
}


print_loss_outcome_BARTLReg <- function(outcome,
                               ct00,ct01,ct10,ct11,
                               cc00,cc01,cc10,cc11,rho){
  d.true <- compute_bivariate_prob_decision(p1_ = outcome[[1]],p_1 = outcome[[2]],
                                            ct00,ct01,ct10,ct11,
                                            cc00,cc01,cc10,cc11,rho)
  d.m1 <- (outcome[[3]] > 0.5)*1
  d.m2 <- (t(outcome[[4]][21,,])>0.5)*1
  
  d.m3 <- compute_bivariate_prob_decision(p1_ = t(outcome[[5]][21,,]),p_1 = t(outcome[[6]][21,,]),
                                          ct00,ct01,ct10,ct11,
                                          cc00,cc01,cc10,cc11,rho)
  
  m1.loss <- d.loss.2(d.m1, 
                      p1_ = outcome[[1]],p_1 = outcome[[2]],rho,
                      ct00,ct01,ct10,ct11,
                      cc00,cc01,cc10,cc11)
  
  m2.loss <- d.loss.2(d.m2, 
                      p1_ = outcome[[1]],p_1 = outcome[[2]],rho,
                      ct00,ct01,ct10,ct11,
                      cc00,cc01,cc10,cc11)
  
  m3.loss <- d.loss.2(d.m3, 
                      p1_ = outcome[[1]],p_1 = outcome[[2]],rho,
                      ct00,ct01,ct10,ct11,
                      cc00,cc01,cc10,cc11)
  
  m1.outcome <- d.outcome.2(d.m1,p1_ = outcome[[1]],p_1 = outcome[[2]])
  m2.outcome <- d.outcome.2(d.m2,p1_ = outcome[[1]],p_1 = outcome[[2]])
  m3.outcome <- d.outcome.2(d.m3,p1_ = outcome[[1]],p_1 = outcome[[2]])
  
  return(list(m1.loss,m2.loss,m3.loss,
              m1.outcome,m2.outcome,m3.outcome))

}

print_acc_metrics_BARTLReg <- function(outcome,
                                        ct00,ct01,ct10,ct11,
                                        cc00,cc01,cc10,cc11,rho){
  d.true <- compute_bivariate_prob_decision(p1_ = outcome[[1]],p_1 = outcome[[2]],
                                            ct00,ct01,ct10,ct11,
                                            cc00,cc01,cc10,cc11,rho)
  d.m1 <- (outcome[[3]] > 0.5)*1
  d.m2 <- (t(outcome[[4]][21,,])>0.5)*1
  
  d.m3 <- compute_bivariate_prob_decision(p1_ = t(outcome[[5]][21,,]),p_1 = t(outcome[[6]][21,,]),
                                          ct00,ct01,ct10,ct11,
                                          cc00,cc01,cc10,cc11,rho)
  
  m1.precision <- d.precision(d.true, d.m1)
  m1.recall <- d.recall(d.true, d.m1)
  m1.f1 <- d.f1(d.true, d.m1)
  
  
  m2.precision <- d.precision(d.true, d.m2)
  m2.recall <- d.recall(d.true, d.m2)
  m2.f1 <- d.f1(d.true, d.m2)
  
  
  m3.precision <- d.precision(d.true, d.m3)
  m3.recall <- d.recall(d.true, d.m3)
  m3.f1 <- d.f1(d.true, d.m3)

  
  return(list(m1.precision,m2.precision,m3.precision,
              m1.recall,m2.recall,m3.recall,
              m1.f1, m2.f1, m3.f1))
}


plot_BART_LReg_loss_outcome <- function(threhold.max, step, filename,scenario,cc00,cc01,cc10,cc11,rho){
  m1.loss <- m2.loss <- m3.loss <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.loss) <- colnames(m2.loss) <- colnames(m3.loss) <- c("weight","mean","lower","upper")
  m1.outcome <- m2.outcome <- m3.outcome <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.outcome) <- colnames(m2.outcome) <- colnames(m3.outcome) <- c("weight","mean","lower","upper")
  i <- 0
  for(threhold in seq(0,threhold.max,step)){
    Filename <- glue::glue(filename,threhold)
    i <- i+1
    ct00 <- 1+threhold
    ct01 <- 1+threhold
    ct10 <- threhold
    ct11 <- threhold
    
    # load the data, which has a list called "outcome"
    load(glue::glue(Filename,scenario,'.RData'))
    loss_outcome <- print_loss_outcome_BARTLReg(outcome,ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
    
    threhold <- threhold*100
    
    m1.loss[i,"weight"] <- threhold
    m1.loss[i,"mean"] <- loss_outcome[[1]][1]
    m1.loss[i,"lower"] <- loss_outcome[[1]][2]
    m1.loss[i,"upper"] <- loss_outcome[[1]][3]
    
    m2.loss[i,"weight"] <- threhold
    m2.loss[i,"mean"] <- loss_outcome[[2]][1]
    m2.loss[i,"lower"] <- loss_outcome[[2]][2]
    m2.loss[i,"upper"] <- loss_outcome[[2]][3]
    
    m3.loss[i,"weight"] <- threhold
    m3.loss[i,"mean"] <- loss_outcome[[3]][1]
    m3.loss[i,"lower"] <- loss_outcome[[3]][2]
    m3.loss[i,"upper"] <- loss_outcome[[3]][3]
    
    m1.outcome[i,"weight"] <- threhold
    m1.outcome[i,"mean"] <- loss_outcome[[4]][1]
    m1.outcome[i,"lower"] <- loss_outcome[[4]][2]
    m1.outcome[i,"upper"] <- loss_outcome[[4]][3]
    
    m2.outcome[i,"weight"] <- threhold
    m2.outcome[i,"mean"] <- loss_outcome[[5]][1]
    m2.outcome[i,"lower"] <- loss_outcome[[5]][2]
    m2.outcome[i,"upper"] <- loss_outcome[[5]][3]
    
    m3.outcome[i,"weight"] <- threhold
    m3.outcome[i,"mean"] <- loss_outcome[[6]][1]
    m3.outcome[i,"lower"] <- loss_outcome[[6]][2]
    m3.outcome[i,"upper"] <- loss_outcome[[6]][3]
    
  }
  
  colors <- c("f1" = "black","f2" = "red", "f3" = "blue")
  shapes <- c("f1" = 1, "f2" = 0, "f3" = 2)

  p.loss <- ggplot() +
    geom_point(aes(x=m1.loss$weight,y=m1.loss$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.loss$weight,y=m2.loss$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.loss$weight,y=m3.loss$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.loss$weight,y=m1.loss$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.loss$weight,y=m2.loss$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.loss$weight,y=m3.loss$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.loss$weight, ymin=m1.loss$lower, ymax=m1.loss$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.loss$weight, ymin=m2.loss$lower, ymax=m2.loss$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.loss$weight, ymin=m3.loss$lower, ymax=m3.loss$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,threhold.max*100,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d^paste("*")), hat(d^paste("*"))[s], hat(d^paste("*"))[s^paste("'")]))+
    labs(x = "decision threshold",
         y = "R(d)",
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=30),
          legend.title=element_text(size=30))
  
  
  p.outcome <- ggplot() +
    geom_point(aes(x=m1.outcome$weight,y=m1.outcome$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.outcome$weight,y=m2.outcome$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.outcome$weight,y=m3.outcome$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.outcome$weight,y=m1.outcome$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.outcome$weight,y=m2.outcome$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.outcome$weight,y=m3.outcome$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.outcome$weight, ymin=m1.outcome$lower, ymax=m1.outcome$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.outcome$weight, ymin=m2.outcome$lower, ymax=m2.outcome$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.outcome$weight, ymin=m3.outcome$lower, ymax=m3.outcome$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,threhold.max*100,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d^paste("*")), hat(d^paste("*"))[s], hat(d^paste("*"))[s^paste("'")]))+
    labs(x = "decision threshold",
         y = "V(d)",
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=30),
          legend.title=element_text(size=30))
  
  
  return(list(p.loss,p.outcome))
}



loss_nsim_nsample <- function(d.m1,L1,L0){
  m1.loss <- apply(L1*d.m1 + L0*(1-d.m1),1,mean)
  m1.loss.mean <- mean(m1.loss)
  m1.loss.sd <- sd(m1.loss)
  m1.loss.upper <- m1.loss.mean + 1.98*m1.loss.sd/sqrt(100)
  m1.loss.lower <- m1.loss.mean - 1.98*m1.loss.sd/sqrt(100)
  m1.loss.sd <- round(m1.loss.sd, digits = 3)
  return(c(m1.loss.mean,m1.loss.lower,m1.loss.upper, m1.loss.sd))
}


plot_BART_TREE_loss_outcome <- function(threhold.max, step, depth, outcome){
  m1.loss <- m2.loss <- m3.loss <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.loss) <- colnames(m2.loss) <- colnames(m3.loss) <- c("weight","mean","lower","upper")
  m1.outcome <- m2.outcome <- m3.outcome <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.outcome) <- colnames(m2.outcome) <- colnames(m3.outcome) <- c("weight","mean","lower","upper")
  m2 <- t(outcome[[5]][depth,,])
  m3 <- t(outcome[[6]][depth,,]-outcome[[7]][depth,,])
  i <- 0
  for(threhold in seq(0,threhold.max,step)){
    i <- i+1

    d.true <- (outcome[[1]]-outcome[[2]] > threhold)*1
    d.m1 <- (outcome[[3]]-outcome[[4]] > threhold)*1 
    d.m2 <- (m2 > threhold)*1
    d.m3 <- (m3 > threhold)*1
    
    p11 <- 0*sqrt(outcome[[1]]*(1-outcome[[1]])*outcome[[2]]*(1-outcome[[2]])) + outcome[[1]]*outcome[[2]]
    p10 <- outcome[[1]] - p11
    p01 <- outcome[[2]] - p11
    p00 <- 1 - p11 - p10 - p01
    L1 <- p00*(1+threhold) + p01*(1+threhold) +p10*threhold + p11*threhold
    L0 <- p00 + p10
    
    threhold <- threhold*100
    m1.loss.value <- loss_nsim_nsample(d.m1,L1,L0)
    m1.loss[i,"weight"] <- threhold
    m1.loss[i,"mean"] <- m1.loss.value[1]
    m1.loss[i,"lower"] <- m1.loss.value[2]
    m1.loss[i,"upper"] <- m1.loss.value[3]
 
    m2.loss.value <- loss_nsim_nsample(d.m2,L1,L0)
    m2.loss[i,"weight"] <- threhold
    m2.loss[i,"mean"] <- m2.loss.value[1]
    m2.loss[i,"lower"] <- m2.loss.value[2]
    m2.loss[i,"upper"] <- m2.loss.value[3]
    
    m3.loss.value <- loss_nsim_nsample(d.m3,L1,L0)
    m3.loss[i,"weight"] <- threhold
    m3.loss[i,"mean"] <- m3.loss.value[1]
    m3.loss[i,"lower"] <- m3.loss.value[2]
    m3.loss[i,"upper"] <- m3.loss.value[3]

    m1.outcome.value <- d.outcome.2(d.m1,outcome[[1]],outcome[[2]])
    m1.outcome[i,"weight"] <- threhold
    m1.outcome[i,"mean"] <- m1.outcome.value[1]
    m1.outcome[i,"lower"] <- m1.outcome.value[2]
    m1.outcome[i,"upper"] <- m1.outcome.value[3]
    
    m2.outcome.value <- d.outcome.2(d.m2,outcome[[1]],outcome[[2]])
    m2.outcome[i,"weight"] <- threhold
    m2.outcome[i,"mean"] <- m2.outcome.value[1]
    m2.outcome[i,"lower"] <- m2.outcome.value[2]
    m2.outcome[i,"upper"] <- m2.outcome.value[3]
    
    m3.outcome.value <- d.outcome.2(d.m3,outcome[[1]],outcome[[2]])
    m3.outcome[i,"weight"] <- threhold
    m3.outcome[i,"mean"] <- m3.outcome.value[1]
    m3.outcome[i,"lower"] <- m3.outcome.value[2]
    m3.outcome[i,"upper"] <- m3.outcome.value[3]
  }
  
  colors <- c("f1" = "black","f2" = "red", "f3" = "blue")
  shapes <- c("f1" = 1, "f2" = 0, "f3" = 2)
  p.loss <- ggplot() +
    geom_point(aes(x=m1.loss$weight,y=m1.loss$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.loss$weight,y=m2.loss$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.loss$weight,y=m3.loss$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.loss$weight,y=m1.loss$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.loss$weight,y=m2.loss$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.loss$weight,y=m3.loss$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.loss$weight, ymin=m1.loss$lower, ymax=m1.loss$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.loss$weight, ymin=m2.loss$lower, ymax=m2.loss$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.loss$weight, ymin=m3.loss$lower, ymax=m3.loss$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,threhold.max*100,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d^paste("*")), hat(d^paste("*"))[s], hat(d^paste("*"))[s^paste("'")]))+
    labs(x = "decision threshold",
         y = "R(d)",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  
  p.outcome <- ggplot() +
    geom_point(aes(x=m1.outcome$weight,y=m1.outcome$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.outcome$weight,y=m2.outcome$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.outcome$weight,y=m3.outcome$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.outcome$weight,y=m1.outcome$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.outcome$weight,y=m2.outcome$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.outcome$weight,y=m3.outcome$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.outcome$weight, ymin=m1.outcome$lower, ymax=m1.outcome$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.outcome$weight, ymin=m2.outcome$lower, ymax=m2.outcome$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.outcome$weight, ymin=m3.outcome$lower, ymax=m3.outcome$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,threhold.max*100,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d^paste("*")), hat(d^paste("*"))[s], hat(d^paste("*"))[s^paste("'")]))+
    labs(x = "decision threshold",
         y = "V(d)",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  return(list(p.loss,p.outcome))
}

plot_BART_TREE_acc_precision_recall <- function(threhold.max, step, depth, outcome){
  #browser()
  m1.acc <- m2.acc <- m3.acc <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.acc) <- colnames(m2.acc) <- colnames(m3.acc) <- c("weight","mean","lower","upper")
  m1.precision <- m2.precision <- m3.precision <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.precision) <- colnames(m2.precision) <- colnames(m3.precision) <- c("weight","mean","lower","upper")
  m1.recall <- m2.recall <- m3.recall <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.recall) <- colnames(m2.recall) <- colnames(m3.recall) <- c("weight","mean","lower","upper")
  m2 <- t(outcome[[5]][depth,,])
  m3 <- t(outcome[[6]][depth,,]-outcome[[7]][depth,,])
  i <- 0
  for(threhold in seq(0,threhold.max,step)){
    i <- i+1

    # Filename <- glue::glue(filename,threhold)
    # i <- i+1
    # ct00 <- 1+threhold
    # ct01 <- 1+threhold
    # ct10 <- threhold
    # ct11 <- threhold
    # 
    # # load the data, which has a list called "outcome"
    # load(glue::glue(Filename,scenario,'.RData'))
    # decisions <- print_d_BARTLReg(outcome,ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
    # d.true <- decisions[[1]]
    # d.m1 <- decisions[[2]]
    
    
    d.true <- ((outcome[[1]]-outcome[[2]]) > threhold)*1
    d.m1 <- ((outcome[[3]]-outcome[[4]]) > threhold)*1 
    d.m2 <- (outcome[[5]] > threhold)*1
    d.m3 <- ((outcome[[6]]-outcome[[7]]) > threhold)*1
    n <- dim(d.true)[1]
    
    threhold <- threhold*100
    
    m1.acc[i,"weight"] <- threhold
    d.m1.acc.mean <- mean(apply((d.m1==d.true)*1, 1, mean))
    d.m1.acc.sd <- sd(apply((d.m1==d.true)*1, 1, mean))
    m1.acc[i,"mean"] <- d.m1.acc.mean
    m1.acc[i,"upper"] <- d.m1.acc.mean + 1.98*d.m1.acc.sd/sqrt(100)
    m1.acc[i,"lower"] <- d.m1.acc.mean - 1.98*d.m1.acc.sd/sqrt(100)
    
    d.m.estime.i <- t(d.m2[depth,,])
    d.m.mean.4sim <- apply((d.m.estime.i == d.true)*1, 1, mean)
    mean <- mean(d.m.mean.4sim)
    se <- sd(d.m.mean.4sim)
    m2.acc[i,"weight"] <- threhold
    m2.acc[i,"mean"] <- mean
    m2.acc[i,"upper"] <- mean + 1.98*se/sqrt(n)
    m2.acc[i,"lower"] <- mean - 1.98*se/sqrt(n)

    d.m.estime.i <- t(d.m3[depth,,])
    d.m.mean.4sim <- apply((d.m.estime.i == d.true)*1, 1, mean)
    mean <- mean(d.m.mean.4sim)
    se <- sd(d.m.mean.4sim)
    m3.acc[i,"weight"] <- threhold
    m3.acc[i,"mean"] <- mean
    m3.acc[i,"upper"] <- mean + 1.98*se/sqrt(n)
    m3.acc[i,"lower"] <- mean - 1.98*se/sqrt(n)

    d.m1.true.positive <- apply(((d.m1==d.true)&(d.true==1))*1,1,sum)
    d.m1.positive <- apply(d.m1,1,sum)
    d.m1.precision <- d.m1.true.positive/d.m1.positive
    d.m1.precision.mean <- mean(na.omit(d.m1.precision))
    d.m1.precision.sd <- sd(na.omit(d.m1.precision))
    m1.precision[i,"weight"] <- threhold
    m1.precision[i,"mean"] <- d.m1.precision.mean
    m1.precision[i,"upper"] <- d.m1.precision.mean+1.98*d.m1.precision.sd/sqrt(100)
    m1.precision[i,"lower"] <- d.m1.precision.mean-1.98*d.m1.precision.sd/sqrt(100)

    d.m.estime.i <- t(d.m2[depth,,])
    d.m.estime.predict <- apply((d.m.estime.i == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.true) & (d.true==1))*1,1,sum)
    precision <- d.m.estime.true/d.m.estime.predict
    mean <- mean(na.omit(precision))
    sd <- sd(na.omit(precision))
    m2.precision[i,"weight"] <- threhold
    m2.precision[i,"mean"] <- mean
    m2.precision[i,"upper"] <- mean + 1.98*sd/sqrt(sum(!is.na(precision)))
    m2.precision[i,"lower"] <- mean - 1.98*sd/sqrt(sum(!is.na(precision)))
    
    d.m.estime.i <- t(d.m3[depth,,])
    d.m.estime.predict <- apply((d.m.estime.i == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.true) & (d.true==1))*1,1,sum)
    precision <- d.m.estime.true/d.m.estime.predict
    mean <- mean(na.omit(precision))
    sd <- sd(na.omit(precision))
    m3.precision[i,"weight"] <- threhold
    m3.precision[i,"mean"] <- mean
    m3.precision[i,"upper"] <- mean + 1.98*sd/sqrt(sum(!is.na(precision)))
    m3.precision[i,"lower"] <- mean - 1.98*sd/sqrt(sum(!is.na(precision)))

  
    d.m1.true.positive <- apply(((d.m1==d.true)&(d.true==1))*1,1,sum)
    d.true.positive <- apply(d.true,1,sum)
    d.m1.recall <- d.m1.true.positive/d.true.positive
    d.m1.recall.mean <- mean(na.omit(d.m1.recall))
    d.m1.recall.sd <- sd(na.omit(d.m1.recall))
    d.m1.recall.upper <- d.m1.recall.mean+1.98*d.m1.recall.sd/sqrt(100)
    d.m1.recall.lower <- d.m1.recall.mean-1.98*d.m1.recall.sd/sqrt(100)
    m1.recall[i,"weight"] <- threhold
    m1.recall[i,"mean"] <- d.m1.recall.mean
    m1.recall[i,"upper"] <- d.m1.recall.upper
    m1.recall[i,"lower"] <- d.m1.recall.lower
    
    
    d.m.estime.i <- t(d.m2[depth,,])
    d.m.true.positive <- apply((d.true == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.true) & (d.true==1))*1,1,sum)
    recall <- d.m.estime.true/d.m.true.positive
    mean <- mean(na.omit(recall))
    sd <- sd(na.omit(recall))
    m2.recall[i,"weight"] <- threhold
    m2.recall[i,"mean"] <- mean
    m2.recall[i,"upper"] <- mean + 1.98*sd/sum(!is.na(recall))
    m2.recall[i,"lower"] <- mean - 1.98*sd/sum(!is.na(recall))
    
    d.m.estime.i <- t(d.m3[depth,,])
    d.m.true.positive <- apply((d.true == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.true) & (d.true==1))*1,1,sum)
    recall <- d.m.estime.true/d.m.true.positive
    mean <- mean(na.omit(recall))
    sd <- sd(na.omit(recall))
    m3.recall[i,"weight"] <- threhold
    m3.recall[i,"mean"] <- mean
    m3.recall[i,"upper"] <- mean + 1.98*sd/sum(!is.na(recall))
    m3.recall[i,"lower"] <- mean - 1.98*sd/sum(!is.na(recall))
  }
  
  colors <- c("f1" = "black","f2" = "red", "f3" = "blue")
  shapes <- c("f1" = 1, "f2" = 0, "f3" = 2)
  p.acc <- ggplot() +
    geom_point(aes(x=m1.acc$weight,y=m1.acc$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.acc$weight,y=m2.acc$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.acc$weight,y=m3.acc$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.acc$weight,y=m1.acc$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.acc$weight,y=m2.acc$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.acc$weight,y=m3.acc$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.acc$weight, ymin=m1.acc$lower, ymax=m1.acc$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.acc$weight, ymin=m2.acc$lower, ymax=m2.acc$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.acc$weight, ymin=m3.acc$lower, ymax=m3.acc$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,threhold.max*100,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d^paste("*")), hat(d^paste("*"))[s], hat(d^paste("*"))[s^paste("'")]))+
    labs(x = "decision threshold",
         y = "accuracy",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  m1.precision <- m1.precision[!is.nan(m1.precision$mean),]
  m2.precision <- m2.precision[!is.nan(m2.precision$mean),]
  m3.precision <- m3.precision[!is.nan(m3.precision$mean),]
  precision.weight.limit <- max(max(m1.precision[,'weight']), max(m2.precision[,'weight']),max(m3.precision[,'weight']))
  p.precision <- ggplot() +
    geom_point(aes(x=m1.precision$weight,y=m1.precision$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.precision$weight,y=m2.precision$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.precision$weight,y=m3.precision$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.precision$weight,y=m1.precision$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.precision$weight,y=m2.precision$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.precision$weight,y=m3.precision$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.precision$weight, ymin=m1.precision$lower, ymax=m1.precision$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.precision$weight, ymin=m2.precision$lower, ymax=m2.precision$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.precision$weight, ymin=m3.precision$lower, ymax=m3.precision$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,precision.weight.limit,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d^paste("*")), hat(d^paste("*"))[s], hat(d^paste("*"))[s^paste("'")]))+
    labs(x = "decision threshold",
         y = "precision",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  m1.recall <- m1.recall[!is.nan(m1.recall$mean),]
  m2.recall <- m2.recall[!is.nan(m2.recall$mean),]
  m3.recall <- m3.recall[!is.nan(m3.recall$mean),]
  recall.weight.limit <- max(max(m1.recall$weight, max(m2.recall$weight),max(m3.recall$weight)))
  if(recall.weight.limit==1) {
    error.bar.width <- 0.5 } else {error.bar.width <- 2}
  p.recall <- ggplot() +
    geom_point(aes(x=m1.recall$weight,y=m1.recall$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.recall$weight,y=m2.recall$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.recall$weight,y=m3.recall$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.recall$weight,y=m1.recall$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.recall$weight,y=m2.recall$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.recall$weight,y=m3.recall$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.recall$weight, ymin=m1.recall$lower, ymax=m1.recall$upper, color="f1"), width=error.bar.width, size=1)+
    geom_errorbar(aes(x=m2.recall$weight, ymin=m2.recall$lower, ymax=m2.recall$upper, color="f2"), width=error.bar.width, size=1)+
    geom_errorbar(aes(x=m3.recall$weight, ymin=m3.recall$lower, ymax=m3.recall$upper, color="f3"), width=error.bar.width, size=1)+
    scale_x_continuous(breaks = seq(0,recall.weight.limit,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d^paste("*")), hat(d^paste("*"))[s], hat(d^paste("*"))[s^paste("'")]))+
    labs(x = "decision threshold",
         y = "recall",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  return(list(p.acc,p.precision,p.recall))
}

get_d_m1_from_BARTLReg_simulation_data <- function(filename,scenario, threhold,cc00,cc01,cc10,cc11,rho){
  # m1 is from BARTLReg simulation data
  Filename <- glue::glue(filename,threhold)
  ct00 <- 1+threhold
  ct01 <- 1+threhold
  ct10 <- threhold
  ct11 <- threhold
  # load the data, which has a list called "outcome"
  load(glue::glue(Filename,scenario,'.RData'))
  decisions <- print_d_BARTLReg(outcome,ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
  d.true <- decisions[[1]]
  d.m1 <- decisions[[2]]
  
  m1.precision <- d.precision(d.true, d.m1)
  m1.recall <- d.recall(d.true, d.m1)
  
  return(list(d.m1, m1.precision, m1.recall))
}

get_d_true_from_BARTLReg_simulation_data <- function(filename,scenario, threhold,cc00,cc01,cc10,cc11,rho){
  # m1 is from BARTLReg simulation data
  Filename <- glue::glue(filename,threhold)
  ct00 <- 1+threhold
  ct01 <- 1+threhold
  ct10 <- threhold
  ct11 <- threhold
  # load the data, which has a list called "outcome"
  load(glue::glue(Filename,scenario,'.RData'))
  decisions <- print_d_BARTLReg(outcome,ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
  d.true <- decisions[[1]]
  return(d.true)
}

plot_BART_TREE_acc_precision_recall_m1_from_BARTLReg_simulation <- function(threhold.max, step, depth, outcome, filename,scenario,ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho){
  #browser()
  m1.acc <- m2.acc <- m3.acc <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.acc) <- colnames(m2.acc) <- colnames(m3.acc) <- c("weight","mean","lower","upper")
  m1.precision <- m2.precision <- m3.precision <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.precision) <- colnames(m2.precision) <- colnames(m3.precision) <- c("weight","mean","lower","upper")
  m1.recall <- m2.recall <- m3.recall <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.recall) <- colnames(m2.recall) <- colnames(m3.recall) <- c("weight","mean","lower","upper")
  m2 <- t(outcome[[5]][depth,,])
  m3 <- t(outcome[[6]][depth,,]-outcome[[7]][depth,,])
  i <- 0
  for(threhold in seq(0,threhold.max,step)){
    i <- i+1
    
    d.true <- get_d_true_from_BARTLReg_simulation_data(filename,scenario,threhold,cc00,cc01,cc10,cc11,rho)
    m1 <- get_d_m1_from_BARTLReg_simulation_data(filename,scenario,threhold,cc00,cc01,cc10,cc11,rho)
    d.m1 <- m1[[1]]
    d.m2 <- (outcome[[5]] > threhold)*1
    d.m3 <- ((outcome[[6]]-outcome[[7]]) > threhold)*1
    n <- dim(d.true)[1]

    threhold <- threhold*100
    
    m1.acc[i,"weight"] <- threhold
    d.m1.acc.mean <- mean(apply((d.m1==d.true)*1, 1, mean))
    d.m1.acc.sd <- sd(apply((d.m1==d.true)*1, 1, mean))
    m1.acc[i,"mean"] <- d.m1.acc.mean
    m1.acc[i,"upper"] <- d.m1.acc.mean + 1.98*d.m1.acc.sd/sqrt(100)
    m1.acc[i,"lower"] <- d.m1.acc.mean - 1.98*d.m1.acc.sd/sqrt(100)
    
    d.m.estime.i <- t(d.m2[depth,,])
    d.m.mean.4sim <- apply((d.m.estime.i == d.true)*1, 1, mean)
    mean <- mean(d.m.mean.4sim)
    se <- sd(d.m.mean.4sim)
    m2.acc[i,"weight"] <- threhold
    m2.acc[i,"mean"] <- mean
    m2.acc[i,"upper"] <- mean + 1.98*se/sqrt(n)
    m2.acc[i,"lower"] <- mean - 1.98*se/sqrt(n)
    
    d.m.estime.i <- t(d.m3[depth,,])
    d.m.mean.4sim <- apply((d.m.estime.i == d.true)*1, 1, mean)
    mean <- mean(d.m.mean.4sim)
    se <- sd(d.m.mean.4sim)
    m3.acc[i,"weight"] <- threhold
    m3.acc[i,"mean"] <- mean
    m3.acc[i,"upper"] <- mean + 1.98*se/sqrt(n)
    m3.acc[i,"lower"] <- mean - 1.98*se/sqrt(n)
    
    m1.precision[i,"weight"] <- threhold
    m1.precision[i,"mean"] <- m1[[2]][1]
    m1.precision[i,"upper"] <-  m1[[2]][2]
    m1.precision[i,"lower"] <- m1[[2]][3]
    
    d.m.estime.i <- t(d.m2[depth,,])
    d.m.estime.predict <- apply((d.m.estime.i == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.true) & (d.true==1))*1,1,sum)
    precision <- d.m.estime.true/d.m.estime.predict
    mean <- mean(na.omit(precision))
    sd <- sd(na.omit(precision))
    m2.precision[i,"weight"] <- threhold
    m2.precision[i,"mean"] <- mean
    m2.precision[i,"upper"] <- mean + 1.98*sd/sqrt(sum(!is.na(precision)))
    m2.precision[i,"lower"] <- mean - 1.98*sd/sqrt(sum(!is.na(precision)))
    
    d.m.estime.i <- t(d.m3[depth,,])
    d.m.estime.predict <- apply((d.m.estime.i == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.true) & (d.true==1))*1,1,sum)
    precision <- d.m.estime.true/d.m.estime.predict
    mean <- mean(na.omit(precision))
    sd <- sd(na.omit(precision))
    m3.precision[i,"weight"] <- threhold
    m3.precision[i,"mean"] <- mean
    m3.precision[i,"upper"] <- mean + 1.98*sd/sqrt(sum(!is.na(precision)))
    m3.precision[i,"lower"] <- mean - 1.98*sd/sqrt(sum(!is.na(precision)))
    
    m1.recall[i,"weight"] <- threhold
    m1.recall[i,"mean"] <- m1[[3]][1]
    m1.recall[i,"upper"] <- m1[[3]][2]
    m1.recall[i,"lower"] <- m1[[3]][3]
    
    d.m.estime.i <- t(d.m2[depth,,])
    d.m.true.positive <- apply((d.true == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.true) & (d.true==1))*1,1,sum)
    recall <- d.m.estime.true/d.m.true.positive
    mean <- mean(na.omit(recall))
    sd <- sd(na.omit(recall))
    m2.recall[i,"weight"] <- threhold
    m2.recall[i,"mean"] <- mean
    m2.recall[i,"upper"] <- mean + 1.98*sd/sum(!is.na(recall))
    m2.recall[i,"lower"] <- mean - 1.98*sd/sum(!is.na(recall))
    
    d.m.estime.i <- t(d.m3[depth,,])
    d.m.true.positive <- apply((d.true == 1)*1,1,sum)
    d.m.estime.true <- apply(((d.m.estime.i == d.true) & (d.true==1))*1,1,sum)
    recall <- d.m.estime.true/d.m.true.positive
    mean <- mean(na.omit(recall))
    sd <- sd(na.omit(recall))
    m3.recall[i,"weight"] <- threhold
    m3.recall[i,"mean"] <- mean
    m3.recall[i,"upper"] <- mean + 1.98*sd/sum(!is.na(recall))
    m3.recall[i,"lower"] <- mean - 1.98*sd/sum(!is.na(recall))
  }
  
  colors <- c("f1" = "black","f2" = "red", "f3" = "blue")
  shapes <- c("f1" = 1, "f2" = 0, "f3" = 2)
  p.acc <- ggplot() +
    geom_point(aes(x=m1.acc$weight,y=m1.acc$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.acc$weight,y=m2.acc$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.acc$weight,y=m3.acc$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.acc$weight,y=m1.acc$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.acc$weight,y=m2.acc$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.acc$weight,y=m3.acc$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.acc$weight, ymin=m1.acc$lower, ymax=m1.acc$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.acc$weight, ymin=m2.acc$lower, ymax=m2.acc$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.acc$weight, ymin=m3.acc$lower, ymax=m3.acc$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,threhold.max*100,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{'sim'}))+
    labs(x = "decision threshold",
         y = "accuracy",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  m1.precision <- m1.precision[!is.nan(m1.precision$mean),]
  m2.precision <- m2.precision[!is.nan(m2.precision$mean),]
  m3.precision <- m3.precision[!is.nan(m3.precision$mean),]
  precision.weight.limit <- max(max(m1.precision[,'weight']), max(m2.precision[,'weight']),max(m3.precision[,'weight']))
  p.precision <- ggplot() +
    geom_point(aes(x=m1.precision$weight,y=m1.precision$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.precision$weight,y=m2.precision$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.precision$weight,y=m3.precision$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.precision$weight,y=m1.precision$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.precision$weight,y=m2.precision$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.precision$weight,y=m3.precision$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.precision$weight, ymin=m1.precision$lower, ymax=m1.precision$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.precision$weight, ymin=m2.precision$lower, ymax=m2.precision$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.precision$weight, ymin=m3.precision$lower, ymax=m3.precision$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,precision.weight.limit,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{'sim'}))+
    labs(x = "decision threshold",
         y = "precision",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  m1.recall <- m1.recall[!is.nan(m1.recall$mean),]
  m2.recall <- m2.recall[!is.nan(m2.recall$mean),]
  m3.recall <- m3.recall[!is.nan(m3.recall$mean),]
  recall.weight.limit <- max(max(m1.recall$weight, max(m2.recall$weight),max(m3.recall$weight)))
  if(recall.weight.limit==1) {
    error.bar.width <- 0.5 } else {error.bar.width <- 2}
  p.recall <- ggplot() +
    geom_point(aes(x=m1.recall$weight,y=m1.recall$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.recall$weight,y=m2.recall$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.recall$weight,y=m3.recall$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.recall$weight,y=m1.recall$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.recall$weight,y=m2.recall$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.recall$weight,y=m3.recall$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.recall$weight, ymin=m1.recall$lower, ymax=m1.recall$upper, color="f1"), width=error.bar.width, size=1)+
    geom_errorbar(aes(x=m2.recall$weight, ymin=m2.recall$lower, ymax=m2.recall$upper, color="f2"), width=error.bar.width, size=1)+
    geom_errorbar(aes(x=m3.recall$weight, ymin=m3.recall$lower, ymax=m3.recall$upper, color="f3"), width=error.bar.width, size=1)+
    scale_x_continuous(breaks = seq(0,recall.weight.limit,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{'sim'}))+
    labs(x = "decision threshold",
         y = "recall",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  return(list(p.acc,p.precision,p.recall))
}


plot_BART_LReg_acc_precision_recall <- function(threhold.max, step,filename,scenario,ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho){
  #browser()
  m1.acc <- m2.acc <- m3.acc <- data.frame(matrix(0,nrow = 0,ncol = 4))
  m1.precision <- m2.precision <- m3.precision <- data.frame(matrix(0,nrow = 0,ncol = 4))
  m1.recall <- m2.recall <- m3.recall <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.acc) <- colnames(m2.acc) <- colnames(m3.acc) <- c("weight","mean","lower","upper")
  colnames(m1.precision) <- colnames(m2.precision) <- colnames(m3.precision) <- c("weight","mean","lower","upper")
  colnames(m1.recall) <- colnames(m2.recall) <- colnames(m3.recall) <- c("weight","mean","lower","upper")
  
  i <- 0
  for(threhold in seq(0,threhold.max,step)){
    Filename <- glue::glue(filename,threhold)
    i <- i+1
    ct00 <- 1+threhold
    ct01 <- 1+threhold
    ct10 <- threhold
    ct11 <- threhold
    
    # load the data, which has a list called "outcome"
    load(glue::glue(Filename,scenario,'.RData'))
    decisions <- print_d_BARTLReg(outcome,ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
    d.true <- decisions[[1]]
    d.m1 <- decisions[[2]]
    d.m2 <- decisions[[3]]
    d.m3 <- decisions[[4]]
    threhold <- threhold*100
    results.scenario <- print_acc_metrics_BARTLReg(outcome,
                                                   ct00,ct01,ct10,ct11,
                                                   cc00,cc01,cc10,cc11,rho)
    
    m1.acc.value <- d.acc(d.m1, d.true)
    m1.acc[i,"weight"] <- threhold
    m1.acc[i,"mean"] <- m1.acc.value[1]
    m1.acc[i,"lower"] <- m1.acc.value[2]
    m1.acc[i,"upper"] <- m1.acc.value[3]
    m1.precision[i,"weight"] <- threhold
    m1.precision[i,"mean"] <- results.scenario[[1]][1]
    m1.precision[i,"upper"] <- results.scenario[[1]][2]
    m1.precision[i,"lower"] <- results.scenario[[1]][3]
    m1.recall[i,"weight"] <- threhold
    m1.recall[i,"mean"] <- results.scenario[[4]][1]
    m1.recall[i,"upper"] <- results.scenario[[4]][2]  
    m1.recall[i,"lower"] <- results.scenario[[4]][3]
    
    m2.acc.value <- d.acc(d.m2, d.true)
    m2.acc[i,"weight"] <- threhold
    m2.acc[i,"mean"] <- m2.acc.value[1]
    m2.acc[i,"lower"] <- m2.acc.value[2]
    m2.acc[i,"upper"] <- m2.acc.value[3]
    m2.precision[i,"weight"] <- threhold
    m2.precision[i,"mean"] <- results.scenario[[2]][1]
    m2.precision[i,"upper"] <- results.scenario[[2]][2]
    m2.precision[i,"lower"] <- results.scenario[[2]][3]
    m2.recall[i,"weight"] <- threhold
    m2.recall[i,"mean"] <- results.scenario[[5]][1]
    m2.recall[i,"upper"] <- results.scenario[[5]][2]  
    m2.recall[i,"lower"] <- results.scenario[[5]][3]
    
    m3.acc.value <- d.acc(d.m3, d.true)
    m3.acc[i,"weight"] <- threhold
    m3.acc[i,"mean"] <- m3.acc.value[1]
    m3.acc[i,"lower"] <- m3.acc.value[2]
    m3.acc[i,"upper"] <- m3.acc.value[3]
    m3.precision[i,"weight"] <- threhold
    m3.precision[i,"mean"] <- results.scenario[[3]][1]
    m3.precision[i,"upper"] <- results.scenario[[3]][2]
    m3.precision[i,"lower"] <- results.scenario[[3]][3]
    m3.recall[i,"weight"] <- threhold
    m3.recall[i,"mean"] <- results.scenario[[6]][1]
    m3.recall[i,"upper"] <- results.scenario[[6]][2]  
    m3.recall[i,"lower"] <- results.scenario[[6]][3]
  }
  colors <- c("f1" = "black","f2" = "red", "f3" = "blue")
  shapes <- c("f1" = 1, "f2" = 0, "f3" = 2)
  p.acc <- ggplot() +
    geom_point(aes(x=m1.acc$weight,y=m1.acc$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.acc$weight,y=m2.acc$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.acc$weight,y=m3.acc$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.acc$weight,y=m1.acc$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.acc$weight,y=m2.acc$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.acc$weight,y=m3.acc$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.acc$weight, ymin=m1.acc$lower, ymax=m1.acc$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.acc$weight, ymin=m2.acc$lower, ymax=m2.acc$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.acc$weight, ymin=m3.acc$lower, ymax=m3.acc$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,threhold.max*100,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d^paste("*")), hat(d^paste("*"))[s], hat(d^paste("*"))[s^paste("'")]))+
    labs(x = "decision threshold",
         y = "accuracy",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  m1.precision <- m1.precision[!is.nan(m1.precision$mean),]
  m2.precision <- m2.precision[!is.nan(m2.precision$mean),]
  m3.precision <- m3.precision[!is.nan(m3.precision$mean),]
  precision.weight.limit <- max(max(m1.precision[,'weight']), max(m2.precision[,'weight']),max(m3.precision[,'weight']))
  p.precision <- ggplot() +
    geom_point(aes(x=m1.precision$weight,y=m1.precision$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.precision$weight,y=m2.precision$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.precision$weight,y=m3.precision$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.precision$weight,y=m1.precision$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.precision$weight,y=m2.precision$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.precision$weight,y=m3.precision$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.precision$weight, ymin=m1.precision$lower, ymax=m1.precision$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.precision$weight, ymin=m2.precision$lower, ymax=m2.precision$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.precision$weight, ymin=m3.precision$lower, ymax=m3.precision$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,precision.weight.limit,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d^paste("*")), hat(d^paste("*"))[s], hat(d^paste("*"))[s^paste("'")]))+
    labs(x = "decision threshold",
         y = "precision",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  m1.recall <- m1.recall[!is.nan(m1.recall$mean),]
  m2.recall <- m2.recall[!is.nan(m2.recall$mean),]
  m3.recall <- m3.recall[!is.nan(m3.recall$mean),]
  recall.weight.limit <- max(max(m1.recall$weight, max(m2.recall$weight),max(m3.recall$weight)))
  p.recall <- ggplot() +
    geom_point(aes(x=m1.recall$weight,y=m1.recall$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.recall$weight,y=m2.recall$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.recall$weight,y=m3.recall$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.recall$weight,y=m1.recall$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.recall$weight,y=m2.recall$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.recall$weight,y=m3.recall$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.recall$weight, ymin=m1.recall$lower, ymax=m1.recall$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.recall$weight, ymin=m2.recall$lower, ymax=m2.recall$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.recall$weight, ymin=m3.recall$lower, ymax=m3.recall$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,recall.weight.limit,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d^paste("*")), hat(d^paste("*"))[s], hat(d^paste("*"))[s^paste("'")]))+
    labs(x = "decision threshold",
         y = "recall",
         color = "Treatment rules",
         shape = "Treatment rules")+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  return(list(p.acc,p.precision,p.recall))
}


d.acc <- function(d.m1, d.true){
  n <- dim(d.m1)[1]
  d.m1.acc.mean <- mean(apply((d.m1==d.true)*1, 1, mean))
  d.m1.acc.sd <- sd(apply((d.m1==d.true)*1, 1, mean))
  d.m1.acc.upper <- d.m1.acc.mean + 1.98*d.m1.acc.sd/sqrt(n)
  d.m1.acc.lower <- d.m1.acc.mean - 1.98*d.m1.acc.sd/sqrt(n)
  return(c(d.m1.acc.mean, d.m1.acc.lower, d.m1.acc.upper))
}


plot_BART_TREE_acc <- function(threhold.max, step, depth, outcome){
  m1.acc <- m2.acc <- m3.acc <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.acc) <- colnames(m2.acc) <- colnames(m3.acc) <- c("weight","mean","lower","upper")
  m2 <- t(outcome[[5]][depth,,])
  m3 <- t(outcome[[6]][depth,,]-outcome[[7]][depth,,])
  i <- 0
  for(threhold in seq(0,threhold.max,step)){
    i <- i+1

    d.true <- (outcome[[1]]-outcome[[2]] > threhold)*1
    d.m1 <- (outcome[[3]]-outcome[[4]] > threhold)*1 
    d.m2 <- (m2 > threhold)*1
    d.m3 <- (m3 > threhold)*1

    threhold <- threhold*100
    
    m1.acc.value <- d.acc(d.m1, d.true)
    m1.acc[i,"weight"] <- threhold
    m1.acc[i,"mean"] <- m1.acc.value[1]
    m1.acc[i,"lower"] <- m1.acc.value[2]
    m1.acc[i,"upper"] <- m1.acc.value[3]
    
    m2.acc.value <- d.acc(d.m2, d.true)
    m2.acc[i,"weight"] <- threhold
    m2.acc[i,"mean"] <- m2.acc.value[1]
    m2.acc[i,"lower"] <- m2.acc.value[2]
    m2.acc[i,"upper"] <- m2.acc.value[3]
    
    m3.acc.value <- d.acc(d.m3, d.true)
    m3.acc[i,"weight"] <- threhold
    m3.acc[i,"mean"] <- m3.acc.value[1]
    m3.acc[i,"lower"] <- m3.acc.value[2]
    m3.acc[i,"upper"] <- m3.acc.value[3]
    
  }
  
  colors <- c("f1" = "black","f2" = "blue", "f3" = "red")
  shapes <- c("f1" = 1, "f2" = 0, "f3" = 2)
  p.acc <- ggplot() +
    geom_point(aes(x=m1.acc$weight,y=m1.acc$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.acc$weight,y=m2.acc$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.acc$weight,y=m3.acc$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.acc$weight,y=m1.acc$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.acc$weight,y=m2.acc$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.acc$weight,y=m3.acc$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.acc$weight, ymin=m1.acc$lower, ymax=m1.acc$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.acc$weight, ymin=m2.acc$lower, ymax=m2.acc$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.acc$weight, ymin=m3.acc$lower, ymax=m3.acc$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,threhold.max*100,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{'sim'}))+
    labs(x = expression(paste(c[t]/c[n],"(%)")),
         y = 'accuracy',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=18),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  return(p.acc)
}


print_d_BARTLReg <- function(outcome,
                               ct00,ct01,ct10,ct11,
                               cc00,cc01,cc10,cc11,rho){
  d.true <- compute_bivariate_prob_decision(p1_ = outcome[[1]],p_1 = outcome[[2]],
                                            ct00,ct01,ct10,ct11,
                                            cc00,cc01,cc10,cc11,rho)
  d.m1 <- (outcome[[3]] > 0.5)*1
  d.m2 <- (t(outcome[[4]][21,,])>0.5)*1
  
  d.m3 <- compute_bivariate_prob_decision(p1_ = t(outcome[[5]][21,,]),p_1 = t(outcome[[6]][21,,]),
                                          ct00,ct01,ct10,ct11,
                                          cc00,cc01,cc10,cc11,rho)

  return(list(d.true,d.m1,d.m2,d.m3))
  
}


plot_BART_LReg_acc <- function(threhold.max, step, filename,scenario,cc00,cc01,cc10,cc11,rho){
  m1.acc <- m2.acc <- m3.acc <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.acc) <- colnames(m2.acc) <- colnames(m3.acc) <- c("weight","mean","lower","upper")
  i <- 0
  for(threhold in seq(0,threhold.max,step)){
    Filename <- glue::glue(filename,threhold)
    i <- i+1
    ct00 <- 1+threhold
    ct01 <- 1+threhold
    ct10 <- threhold
    ct11 <- threhold
    
    # load the data, which has a list called "outcome"
    load(glue::glue(Filename,scenario,'.RData'))
    decisions <- print_d_BARTLReg(outcome,ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
    d.true <- decisions[[1]]
    d.m1 <- decisions[[2]]
    d.m2 <- decisions[[3]]
    d.m3 <- decisions[[4]]
    threhold <- threhold*100
    
    m1.acc.value <- d.acc(d.m1, d.true)
    m1.acc[i,"weight"] <- threhold
    m1.acc[i,"mean"] <- m1.acc.value[1]
    m1.acc[i,"lower"] <- m1.acc.value[2]
    m1.acc[i,"upper"] <- m1.acc.value[3]
    
    m2.acc.value <- d.acc(d.m2, d.true)
    m2.acc[i,"weight"] <- threhold
    m2.acc[i,"mean"] <- m2.acc.value[1]
    m2.acc[i,"lower"] <- m2.acc.value[2]
    m2.acc[i,"upper"] <- m2.acc.value[3]
    
    m3.acc.value <- d.acc(d.m3, d.true)
    m3.acc[i,"weight"] <- threhold
    m3.acc[i,"mean"] <- m3.acc.value[1]
    m3.acc[i,"lower"] <- m3.acc.value[2]
    m3.acc[i,"upper"] <- m3.acc.value[3]
  }
  
  colors <- c("f1" = "black","f2" = "blue", "f3" = "red")
  shapes <- c("f1" = 1, "f2" = 0, "f3" = 2)
  
  p.acc <- ggplot() +
    geom_point(aes(x=m1.acc$weight,y=m1.acc$mean,color="f1"),size=3)+
    geom_point(aes(x=m2.acc$weight,y=m2.acc$mean,color="f2"),size=3)+
    geom_point(aes(x=m3.acc$weight,y=m3.acc$mean,color="f3"),size=3)+
    geom_line(aes(x=m1.acc$weight,y=m1.acc$mean,color="f1"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m2.acc$weight,y=m2.acc$mean,color="f2"),size=0.5,linetype = "dashed")+
    geom_line(aes(x=m3.acc$weight,y=m3.acc$mean,color="f3"),size=0.5,linetype = "dashed")+
    geom_errorbar(aes(x=m1.acc$weight, ymin=m1.acc$lower, ymax=m1.acc$upper, color="f1"), width=2, size=1)+
    geom_errorbar(aes(x=m2.acc$weight, ymin=m2.acc$lower, ymax=m2.acc$upper, color="f2"), width=2, size=1)+
    geom_errorbar(aes(x=m3.acc$weight, ymin=m3.acc$lower, ymax=m3.acc$upper, color="f3"), width=2, size=1)+
    scale_x_continuous(breaks = seq(0,threhold.max*100,step*100))+
    scale_color_manual(values = colors,labels = expression(hat(d)^paste("*"), hat(d)^sim, hat(d)[com]^{'sim'}))+
    labs(x = expression(paste(c[t]/c[n],"(%)")),
         y = 'accuracy',
         color = "Model",
         shape = "Model")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=18),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18))
  
  return(p.acc)
}

get_acc <- function(outcome,threhold,cc00,cc01,cc10,cc11,rho){
  #browser()
  m1.acc <- m2.acc <- m3.acc <- data.frame(matrix(0,nrow = 0,ncol = 4))
  colnames(m1.acc) <- colnames(m2.acc) <- colnames(m3.acc) <- c("weight","mean","lower","upper")
  ct00 <- 1+threhold
  ct01 <- 1+threhold
  ct10 <- threhold
  ct11 <- threhold
  
  i <- 1
  decisions <- print_d_BARTLReg(outcome,ct00,ct01,ct10,ct11,cc00,cc01,cc10,cc11,rho)
  d.true <- decisions[[1]]
  d.m1 <- decisions[[2]]
  d.m2 <- decisions[[3]]
  d.m3 <- decisions[[4]]
  threhold <- threhold*100
    
  m1.acc.value <- d.acc(d.m1, d.true)
  m1.acc[i,"weight"] <- threhold
  m1.acc[i,"mean"] <- m1.acc.value[1]
  m1.acc[i,"lower"] <- m1.acc.value[2]
  m1.acc[i,"upper"] <- m1.acc.value[3]
    
  m2.acc.value <- d.acc(d.m2, d.true)
  m2.acc[i,"weight"] <- threhold
  m2.acc[i,"mean"] <- m2.acc.value[1]
  m2.acc[i,"lower"] <- m2.acc.value[2]
  m2.acc[i,"upper"] <- m2.acc.value[3]
    
  m3.acc.value <- d.acc(d.m3, d.true)
  m3.acc[i,"weight"] <- threhold
  m3.acc[i,"mean"] <- m3.acc.value[1]
  m3.acc[i,"lower"] <- m3.acc.value[2]
  m3.acc[i,"upper"] <- m3.acc.value[3]

  return(list(m1.acc,m2.acc,m3.acc))
}


get_BART_TREE_loss_outcome <- function(threhold, depth, outcome,scenario){
  m1.loss <- m2.loss <- m3.loss <- data.frame(matrix(0,nrow = 0,ncol = 3))
  colnames(m1.loss) <- colnames(m2.loss) <- colnames(m3.loss) <- c("weight","mean","sd")
  m1.outcome <- m2.outcome <- m3.outcome <- data.frame(matrix(0,nrow = 0,ncol = 3))
  colnames(m1.outcome) <- colnames(m2.outcome) <- colnames(m3.outcome) <- c("weight","mean","sd")
  m2 <- t(outcome[[5]][depth,,])
  m3 <- t(outcome[[6]][depth,,]-outcome[[7]][depth,,])
  i <- 1
  
  d.true <- (outcome[[1]]-outcome[[2]] > threhold)*1
  d.m1 <- (outcome[[3]]-outcome[[4]] > threhold)*1 
  d.m2 <- (m2 > threhold) * 1
  d.m3 <- (m3 > threhold) * 1
  
  p11 <- 0 * sqrt(outcome[[1]] * (1 - outcome[[1]]) * outcome[[2]] * (1 - outcome[[2]])) + outcome[[1]] * outcome[[2]]
  p10 <- outcome[[1]] - p11
  p01 <- outcome[[2]] - p11
  p00 <- 1 - p11 - p10 - p01
  L1 <- p00 * (1 + threhold) + p01 * (1 + threhold) + p10 * threhold + p11 * threhold
  L0 <- p00 + p10
  
  threhold <- threhold * 100
  m1.loss.value <- loss_nsim_nsample(d.m1, L1, L0)
  m1.loss[i, "weight"] <- threhold
  m1.loss[i, "mean"] <- m1.loss.value[1]
  m1.loss[i, "sd"] <- m1.loss.value[4]

  m2.loss.value <- loss_nsim_nsample(d.m2, L1, L0)
  m2.loss[i, "weight"] <- threhold
  m2.loss[i, "mean"] <- m2.loss.value[1]
  m2.loss[i, "sd"] <- m2.loss.value[4]

  m3.loss.value <- loss_nsim_nsample(d.m3, L1, L0)
  m3.loss[i, "weight"] <- threhold
  m3.loss[i, "mean"] <- m3.loss.value[1]
  m3.loss[i, "sd"] <- m3.loss.value[4]

  m1.outcome.value <- d.outcome.2(d.m1, outcome[[1]], outcome[[2]])
  m1.outcome[i, "weight"] <- threhold
  m1.outcome[i, "mean"] <- m1.outcome.value[1]
  m1.outcome[i, "sd"] <- m1.outcome.value[4]

  m2.outcome.value <- d.outcome.2(d.m2, outcome[[1]], outcome[[2]])
  m2.outcome[i, "weight"] <- threhold
  m2.outcome[i, "mean"] <- m2.outcome.value[1]
  m2.outcome[i, "sd"] <- m2.outcome.value[4]

  m3.outcome.value <- d.outcome.2(d.m3, outcome[[1]], outcome[[2]])
  m3.outcome[i, "weight"] <- threhold
  m3.outcome[i, "mean"] <- m3.outcome.value[1]
  m3.outcome[i, "sd"] <- m3.outcome.value[4]
  
  m1.loss <- round(m1.loss,digits = 3)
  m2.loss <- round(m2.loss,digits = 3)
  m3.loss <- round(m3.loss,digits = 3)
  
  m1.outcome <- round(m1.outcome,digits = 3)
  m2.outcome <- round(m2.outcome,digits = 3)
  m3.outcome <- round(m3.outcome,digits = 3)
  
  loss <- cbind(scenario=rep(scenario,3),model=c("f1","f2","f3"),rbind(m1.loss, m2.loss, m3.loss))
  outcome <- cbind(scenario=rep(scenario,3),model=c("f1","f2","f3"),rbind(m1.outcome,m2.outcome,m3.outcome))
  
  
  return(list(loss=loss,
              outcome=outcome))
}


# Perform leave one out prediction using single or multiple regression
LeaveOneOutLinearPrediction <- function(df, interact = FALSE, test_reduced = TRUE) {
  # df            A data frame with at least 2 columns, the first column is the dependent variable, 2 to N are independent variable
  # interact      Whether the linear model should include interaction of independent variables
  # test_reduced  Whether to test the contribution of each independent variable with a reduced model
  
  require(stats); 
  require(boot);
  require(MASS); 
  
  for (i in 2:ncol(df)) if (is.character(df[[i]])) df[[i]] <- factor(df[[i]]); 
  
  out <- list(data=df);
  
  cnm <- colnames(df); 
  fm0 <- paste(cnm[1], '~', paste(cnm[-1], collapse=c(' + ', ' * ')[as.integer(interact)+1]));
  fm0 <- as.formula(fm0); 
  
  # Model with all predictors and observations
  fll <- lm(fm0, data=df);
  smm <- summary(fll); 
  fst <- smm$fstatistic;
  stt <- c(smm$r.squared, smm$adj.r.squared, pf(fst[1], fst[2], fst[3], lower.tail = FALSE), extractAIC(fll)[2]);
  ftt <- fitted(fll); 
  ftt <- ftt[rownames(df)];
  names(stt) <- c('r_squared', 'r_squared_adjusted', 'p_value', 'aic');
  names(ftt) <- rownames(df); 
  
  out$full <- list();
  out$full$stats <- stt;
  out$full$model <- fll;
  
  ################################################################################################
  # Leave one out prediction
  run.leave.one.out <- function(df, interact) {
    l1o <- sapply(1:nrow(df), function(i) {
      cnm <- colnames(df); 
      fm0 <- paste(cnm[1], '~', paste(cnm[-1], collapse=c(' + ', ' * ')[as.integer(interact)+1]));
      df1 <- df[-i, ]; 
      lm1 <- lm(fm0, data=df1); 
      prd <- predict(lm1, df[i, ]); 
      smm <- summary(lm1); 
      fst <- smm$fstatistic;
      stt <- c(predicted = as.numeric(prd), 'r.squared'=smm$r.squared, 'adj.r.squared'=smm$adj.r.squared, 
               'p.value'=pf(fst[1], fst[2], fst[3], lower.tail = FALSE), AIC = extractAIC(lm1)[2]);
      names(stt) <- c('predicted', 'r_squared', 'r_squared_adjusted', 'p_value', 'aic');
      stt;
    }); 
    colnames(l1o) <- rownames(df); 
    t(l1o); 
  }
  
  l1o <- run.leave.one.out(df, interact);
  prd <- cbind(observed=df[, 1], predicted=ftt, leave1out = l1o[, 1]); 
  
  # lmg <- glm(fm0, data=df);
  # err <- cv.glm(df[names(lmg$residuals), ], lmg);
  out$leave1out <- list(MSE=mean((prd[,3]-prd[,1])^2, na.rm=TRUE), 
                        N=nrow(prd[!is.na(prd[,1])&!is.na(prd[, 3]), , drop=FALSE]), 
                        prediction=l1o);

  ################################################################################################
  # Reduce model to evaluate contribution of each predictor
  if (test_reduced & ncol(df) > 2) {
    mdl <- lapply(2:ncol(df), function(i) {
      # Create reduced model
      vrb <- cnm[-c(1, i)]; 
      fm1 <- paste(cnm[1], '~', paste(vrb, collapse=c(' + ', ' * ')[as.integer(interact)+1]));
      fm1 <- as.formula(fm1);
      rdc <- lm(fm1, data=df); 
      
      # Compare reduced and full models
      df1 <- df[!is.na(df[, i])&df[,i]!='', ]; 
      lm0 <- lm(fm0, data=df1);
      lm1 <- lm(fm1, data=df1);
      cmp <- anova(lm0, lm1, test='Chisq');

      # summary stats
      aic <- c(extractAIC(lm0)[2], extractAIC(lm1)[2]); 
      pvl <- cmp$`Pr(>Chi)`[2];
      rdf <- cmp$Res.Df[2];
      
      list(model = rdc, reduce = cmp, stat=c(rdf, pvl, aic)); 
    });

    out$reduced <- list();
    out$reduced$stats <- t(sapply(mdl, function(m) m[[3]]));
    out$reduced$model <- lapply(mdl, function(x) x[[1]]); 
    out$reduced$full_vs_reduced <- lapply(mdl, function(x) x[[2]]);
    
    names(out$reduced$model) <- names(out$reduced$full_vs_reduced) <- rownames(out$reduced$stats) <- cnm[-1];
    colnames(out$reduced$stats) <- c('Residual_Df', 'P_Chisq', 'AIC_Full', 'AIC_Reduced');
    
    # leave one out
    l1o <- sapply(2:ncol(df), function(i) run.leave.one.out(df[, -i], interact)[, 1]); 
    colnames(l1o) <- paste('reduced', cnm[-1], sep='_');
    prd <- cbind(prd, l1o); 
  }
  
  out$prediction <- prd; 
  out;
}

# Summarize a bunch of models from the <LeaveOneOutLinearPrediction> function
SummarizeLeaveOneOutLinearPrediction <- function(models) {
  # models  a list of outputs from runs of <LeaveOneOutLinearPrediction> function into a stat table
  
  mdl <- models;
  
  prd <- sapply(mdl, function(x) x$prediction[, 'leave1out']); 
  obs <- sapply(mdl, function(x) x$prediction[, 'observed']); 
  dff <- prd - obs; 
  mer <- as.vector(colMeans(abs(dff), na.rm=TRUE));
  crr <- sapply(1:length(mdl), function(i) {
    r <- cor.test(obs[, i], prd[, i], use='pair');
    c(r$estimate, r$p.value);
  });
  crr <- t(crr);
  mse <- as.vector(sapply(mdl, function(x) x$leave1out$MSE));
  num <- as.vector(sapply(mdl, function(x) x$leave1out$N));
  pks <- as.vector(sapply(1:length(mdl), function(i) ks.test(prd[, i], obs[, i])$p.value[[1]]));
  pfl <- as.vector(sapply(mdl, function(x) x$full$stats['p_value'])); 
  prs <- as.vector(sapply(1:length(mdl), function(i) shapiro.test(prd[, i]-obs[, i])$p.value));
  prc <- as.vector(sapply(mdl, function(x) rev(x$reduced$stats[, 'P_Chisq'])[[1]]));
  stt <- cbind('N'=num, 'MSE'=mse, 'P_Residual'=prs, 'P_KS'=pks,
               'P_Reduced'=prc, 'Corr_LOOCV'=as.vector(crr[, 1]), 'P_Corr_LOOCV'=as.vector(crr[, 2]), 'P_Full'=pfl);
  rownames(stt) <- names(mdl); 
  
  stt; 
}
# Use a stepwise procedure to add predictors one by one until no more improvement of the model
StepForwardLinearPrediction <- function(df0, df1, top=1, p.cutoff=0.05, max.step=10) {
  # df0       A data frame with at least 1 column, the 2 to n columns are confounders 
  # df1       The data matrix of predictors to be selected
  # top       Number of top predictors to select at each step
  # p.cutoff  The cutoff of p values to stop the procedure
  # max.step  Maximum number of predictors allowed
  
  require(stats); 
  require(boot);
  require(MASS); 
  
  top <- max(1, top);
  p.cutoff <- min(max(p.cutoff, 0), 0.25); 
  
  for (i in 2:ncol(df)) if (is.character(df[[i]])) df[[i]] <- factor(df[[i]]); 
  
  npr <- 0; # current number of predictors in the model
  fml <- paste(colnames(df0)[1], '~'); 
  if (ncol(df0) > 1) fml <- paste(fml, paste(colnames(df0)[-1], collapse=' + ')); 
  
  p.min <- -1; # minimal p value of all models
  step  <- 1; 
  survivors  <- list();
  selected   <- 0; 
  
  ###############################################################################################
  # Add the first set of predictors
  fml <- paste(fml, ' + Predictor1', sep='');
  mdl <- lapply(colnames(df1), function(cnm) {
    d <- cbind(df0, Predictor1=df1[, cnm]); 
    lm(fml, data = d); 
  });
  pvl <- sapply(mdl, function(x) summary(x)$coefficients['Predictor1', 4]);
  names(mdl) <- names(pvl) <- colnames(df1);
  mdl <- mdl[pvl<=p.cutoff];
  mdl <- mdl[order(pvl[names(mdl)])];
  if (length(mdl)>0) {
    mdl <- mdl[1:min(top, length(mdl))]; 
    survivors$step1 <- lapply(1:length(mdl), function(i) list(predictor=names(mdl)[i], model=mdl[[i]]));
    names(survivors$step1) <- names(mdl); 
    selected  <- length(mdl); 
  }; 
  
  ###############################################################################################
  # Add predictors step by step
  while (selected > 0 & step <= max.step) {
    step <- step+1;
    cat('Step', step, '\n'); 
    last.model <- survivors[[length(survivors)]];
    pred.name  <- paste('Predictor', step, sep=''); 
    fml        <- paste(fml, ' + Predictor', step, sep='');
    
    mdls <- lapply(1:length(last.model), function(i) {
      prd0 <- last.model[[i]][[1]];
      mdl0 <- last.model[[i]][[2]];
      df2  <- df1[, !(colnames(df1) %in% prd0), drop=FALSE];
      mdl <- lapply(colnames(df2), function(cnm) {
        d <- cbind(mdl0$model, df2[, cnm]); 
        colnames(d)[ncol(d)] <- pred.name; 
        lm(fml, data = d); 
      });
      names(mdl) <- colnames(df2);
      mdl;
    }); 
    pvls <- lapply(mdls, function(mdl) sapply(mdl, function(x) summary(x)$coefficients[pred.name, 4]));
    info <- data.frame(id = unlist(lapply(pvls, names), use.names=FALSE), 
                       initial=rep(1:length(pvls), sapply(pvls, length)), 
                       pvalue=unlist(pvls, use.names=FALSE), stringsAsFactors = FALSE);
    info <- info[order(info[, 3]), , drop=FALSE]; 
    info <- info[info[, 3] <= p.cutoff, , drop=FALSE]; 
    
    if (nrow(info) > 0) {
      info <- info[1:min(top, nrow(info)), , drop=FALSE];
      ind  <- length(survivors) + 1; 
      survivors[[ind]] <- lapply(1:nrow(info), function(i) {
        m <- mdls[[info[i, 2]]][[info[i, 1]]];
        p <- c(last.model[[info[i, 2]]]$predictor, info[i, 1]);
        list(predictor=p, model=m); 
      })
      names(survivors)[ind] <- paste('step', ind, sep=''); 
      names(survivors[[ind]]) <- sapply(survivors[[ind]], function(x) paste(x[[1]], collapse=' + ')); 
    }; 
    
    selected <- nrow(info); 
  }; 

  survivors;  
}
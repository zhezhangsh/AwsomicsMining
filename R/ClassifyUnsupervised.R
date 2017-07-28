# This function run a number of algorithms for the unsupervised clustering of samples into groups.

ClassifyUnsupervised <- function(predictor, num_group=2, 
                                 algorithms=c('hierarchical', 'kmeans', 'pam', 'som', 'supersom', 
                                              'bic', 'icl', 'svc')) {
  # predictor   a data.frame of one predictor per column, used to classify rows
  # num_group   the number of groups for the rows to be classified into
  # algorithms  one or multiple algorithm to be used to perform the classification
    # hierarchical  Hierarchical clustering
    # kmeans        K-means clustering
    # pam           Partitioning Around Medoids clustering
    # bic           Bayesian Information Criterion for Gaussian mixture models
    # icl           Integrated Complete-data Likelihood for Gaussian mixture models 
    # som           Self-organizing map
    # supersom      An extension of self-organizing map
    # svc           Support vector clustering
  
  require(e1071); 
  require(mclust); 
  require(som);
  require(kohonen); 
  require(neuralnet); 
  require(cluster);
  
  if (class(predictor) != 'data.frame') predictor <- as.data.frame(predictor); 
  if (num_group < 2) num_group <- 2 else num_group <- round(num_group); 
  algorithms <- tolower(algorithms);

  typ <- sapply(1:ncol(predictor), function(i) is.numeric(predictor[[i]])); 
  
  out <- list(classified=NULL, algorithm=NULL); # 1st element is a matrix of group assignment; 2nd is output of each algorithm
  res <- list(); 
  
  #################################################################################
  # Hierarchical clustering
  if ('hierarchical' %in% algorithms & length(typ[!typ])==0) {
    hcl <- hclust(as.dist(1-cor(t(as.matrix(predictor))))); 
    prd <- cutree(hcl, k=num_group); 
    res$Hierarchical <- list(prediction=prd, model=hcl); 
  };
  #################################################################################
  # K-means clustering
  if ('kmeans' %in% algorithms & length(typ[!typ])==0) {
    kmn <- kmeans(as.matrix(predictor), num_group); 
    prd <- kmn$cluster;
    res$Kmeans <- list(prediction=prd, model=kmn); 
  };
  #################################################################################
  # PAM (Partitioning Around Medoids) clustering
  if ('pam' %in% algorithms & length(typ[!typ])==0) {
    pam <- pam(predictor, k=num_group);
    prd <- pam$clustering; 
    res$PAM <- list(prediction=prd, model=pam); 
  };
  #################################################################################
  # Bayesian Information Criterion for Gaussian mixture models
  if ('bic' %in% algorithms & length(typ[!typ])==0) {
    bic <- mclustBIC(predictor, G=num_group);
    mnm <- rownames(as.matrix(summary(bic)))[1];
    mnm <- sub(',[0-9]$', '', mnm);
    mdl <- Mclust(predictor, G=num_group, modelNames = mnm); 
    prd <- mdl$classification; 
    res$BIC <- list(prediction=prd, model=mdl); 
  };
  #################################################################################
  # Integrated Complete-data Likelihood for Gaussian mixture models 
  if ('icl' %in% algorithms & length(typ[!typ])==0) {
    icl <- mclustICL(predictor, G=num_group);
    mnm <- rownames(as.matrix(summary(icl)))[1];
    mnm <- sub(',[0-9]$', '', mnm);
    mdl <- Mclust(predictor, G=num_group, modelNames = mnm); 
    prd <- mdl$classification; 
    res$ICL <- list(prediction=prd, model=mdl); 
  };
  #################################################################################
  # Self organizing map
  if ('som' %in% algorithms & length(typ[!typ])==0) {
    mld <- som::som(as.matrix(predictor), 1, num_group); 
    prd <- mld$visual[, 2]+1;
    res$SOM <- list(prediction=prd, model=mdl); 
  }; 
  #################################################################################
  # Super Self organizing map
  if ('supersom' %in% algorithms & length(typ[!typ])==0) {
    mld <- supersom(as.matrix(predictor), grid=somgrid(xdim = 1, ydim = num_group)); 
    prd <- mld$unit.classif;
    res$SuperSOM <- list(prediction=prd, model=mdl); 
  }; 
  #################################################################################
  # Support vector clustering
  # if ('svc' %in% algorithms & num_group==2) {
  #   mdl <- svm(~., data=predictor);
  #   prd <- as.integer(predict(mdl)) + 1;
  #   res$SVC <- list(prediction=prd, model=mdl);
  # };
  #################################################################################  
  
  cls <- sapply(res, function(x) x[[1]]); 
  dimnames(cls) <- list(rownames(predictor), names(res));
}
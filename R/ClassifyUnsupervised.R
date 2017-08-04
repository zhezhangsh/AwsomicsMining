# This function run a number of algorithms for the unsupervised clustering of samples into groups.

ClassifyUnsupervised <- function(
  predictor, num_group=2, 
  algorithms=c('hierarchical', 'kmeans', 'cmeans', 'cshell', 'pam', 'bagged', 'som', 'supersom', 
               'bic', 'icl', 'svc', 'randomforest')) {
  # predictor   a data.frame of one predictor per column, used to classify rows
  # num_group   the number of groups for the rows to be classified into
  # algorithms  one or multiple algorithm to be used to perform the classification
    # hierarchical  Hierarchical clustering
    # kmeans        K-means clustering
    # cmeans        Fuzzy version of the kmeans
    # cshell        Fuzzy C-Shell Clustering
    # pam           Partitioning Around Medoids clustering
    # bagged        Bagged Clustering
    # bic           Bayesian Information Criterion for Gaussian mixture models
    # icl           Integrated Complete-data Likelihood for Gaussian mixture models 
    # som           Self-organizing map
    # supersom      An extension of self-organizing map
    # svc           Support vector clustering
    # randomforest  Random forest to get distance matrix (times 2 samples in the same terminal nodes), followed by kmeans
  
  require(e1071); 
  require(mclust); 
  require(som);
  require(kohonen); 
  require(neuralnet); 
  require(cluster);
  require(randomForest); 
  
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
  # C-means clustering
  if ('cmeans' %in% algorithms & length(typ[!typ])==0) {
    cmn <- cmeans(as.matrix(predictor), num_group); 
    prd <- cmn$cluster;
    res$Cmeans <- list(prediction=prd, model=cmn); 
  };
  #################################################################################
  # Fuzzy C-Shell Clustering
  if ('cshell' %in% algorithms & length(typ[!typ])==0) {
    csh <- cshell(as.matrix(predictor), num_group); 
    prd <- cmn$cluster;
    res$Cshell <- list(prediction=prd, model=csh); 
  };
  #################################################################################
  # Bagged Clustering
  if ('bagged' %in% algorithms & length(typ[!typ])==0) {
    bcl <- bclust(as.matrix(predictor), centers=num_group); 
    prd <- bcl$cluster;
    res$Bagged <- list(prediction=prd, model=bcl); 
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
  if ('svc' %in% algorithms & num_group==2) {
    mdl <- svm(~., data=predictor);
    prd <- as.integer(predict(mdl)) + 1;
    res$SVC <- list(prediction=prd, model=mdl);
  };
  #################################################################################
  # Random forest
  if ('randomforest' %in% algorithms) {
    rft <- randomForest(predictor);
    prd <- kmeans(rft$proximity, num_group)$cluster;
    res$RandomForest <- list(prediction=prd, model=rft); 
  }
  
  cls <- sapply(res, function(x) x[[1]]); 
  dimnames(cls) <- list(rownames(predictor), names(res));
  
  invisible(list(classification=cls, model=lapply(res, function(x)x[[2]])))
}
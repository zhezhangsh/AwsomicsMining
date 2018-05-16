# Identify the nodes of a network that have more than expected connections to a group of nodes of interest
FindKevinBacon <- function(prs, noi, uni=NULL) {
  # prs: a matrix or data.frame of 2 columns, each row is a pair of nodes
  # noi: nodes of interest
  # uni: the universe; if not empty, limit pairs to those with both ends in the universe
  
  if (length(uni) > 0) prs <- prs[prs[, 1] %in% uni & prs[, 2] %in% uni, , drop=FALSE];
    
  uni <- unique(c(prs[, 1], prs[, 2]));
  noi <- noi[noi %in% c(prs[, 1], prs[, 2])]; 
  neg <- setdiff(uni, noi);
  
  mp1 <- split(prs[, 1], prs[, 2]);
  mp2 <- split(prs[, 2], prs[, 1]);
  map <- split(c(mp1, mp2), c(names(mp1), names(mp2)));
  map <- lapply(map, function(m) unique(unlist(m, use.names = FALSE)));
  frq <- 1/sapply(map, length);
  
  w0 <- n0 <- rep(0, length(uni));
  names(w0) <- names(n0) <- uni;
  n0[noi] <- 1;
  
  stat <- sapply(1:length(map), function(i) { 
    if (i/100==round(i/100)) print(i); 
    nm <- names(map)[i];
    mp <- map[[nm]];
    fq <- frq[nm]*frq[mp];
    w1 <- w0;
    w1[mp] <- fq;
    w1 <- split(w1, n0);
    w1 <- w1[c('0', '1')];
    tt <- wilcox.test(w1[[1]], w1[[2]]);
    mn <- sapply(w1, mean);
    c(length(mp), length(intersect(mp, noi)), mn, tt$p.value);
  });
  
  stat <- t(stat);
  rownames(stat) <- names(map);
  colnames(stat) <- c('Total', 'Within', 'Mean_Out', 'Mean_In', 'P_RST');
  stat <- cbind(stat, FDR=p.adjust(stat[, 5], method='BH'));
  stat <- cbind(stat, Enrichment=(stat[,2]/stat[,1])/(length(noi)/length(uni)));
  stat <- stat[, c(1, 2, 7, 3:6)];
  
  stat;
}
require("yaGST")
require("parallel")

#Functions

#============================================================================================
loading_data <- function(filename,M){
  z <- readRDS(file=paste0("Data/",filename,"/TCGA-",filename,"_normcounts.rda"))
  
  #Get the gene names and associate it with appropriate variables
  #====================================================================================
  load("Data/Others/TF_Target_Info.Rdata")
  tfs <- tf_gene_names;
  targets <- target_gene_names;
  rownames(z) <- target_gene_names;

  #Get TFs and targets
  all_genes <- rownames(z);
  m_targets <- colnames(M);
  m_tfs <- rownames(M);
  temp_target_genes <- intersect(all_genes,m_targets);
  tf_genes <- intersect(all_genes,m_tfs);
  rem_target_genes <- setdiff(temp_target_genes,tf_genes);
  target_genes <- c(tf_genes,rem_target_genes);
  
  #Get the modified gene expression matrix
  genes_from_expression_matrix <- union(target_genes,tf_genes);
  modified_exp_matrix <- t(z[genes_from_expression_matrix,]);
  rm(z);
  gc();
  N <- nrow(modified_exp_matrix);
  d <- ncol(modified_exp_matrix);
  
  #Get final mechanistic network
  g_M <- M[tf_genes,target_genes];
  
  #Make the perturbation matrix
  K <- matrix(0,nrow=N,ncol=d);
  colnames(K) <- colnames(modified_exp_matrix);
  
  return(list(modified_exp_matrix,K,g_M,tf_genes,target_genes));
}


#===========================================================================================
parcor <- function(Mat) { 
  
  require(parallel) 
  nc <- detectCores() 
  m <- dim(Mat)[1] 
  n <- dim(Mat)[2] 
  rL <- split(t(Mat),1:n) 
  res <- mclapply(rL, function(x){ mm <- matrix(unlist(x), ncol = m, byrow = TRUE) 
                  cor(t(mm), Mat)
                  },
                  mc.cores = nc ) 
  out <- matrix(unlist(res), nrow=n, byrow = T) 
  return(out)
}

#===========================================================================================
#Get regulons for each candidate master regulator
constructRegulon2 <- function(net, corNet, minReg = 20){
  net <- net[rowSums(net != 0) > minReg, ]
  netReg2 <- vector("list", nrow(net))
  names(netReg2) <- rownames(net)
  toKeep <- NULL
  for(k in 1:nrow(net)){
    tf <- rownames(net)[k]
    #print(tf)
    reg <- colnames(net[tf, net[tf, ] != 0, drop = F])
    reg <- reg[reg %in% colnames(corNet)]
    tmp <- corNet[tf, reg, drop = F]
    
    pos <- colnames(tmp[1, tmp[1, ] > 0, drop = F])
    neg <- colnames(tmp[1, tmp[1, ] < 0, drop = F])
    
    if(length(pos) < minReg & length(neg) < minReg) next
    if(length(pos) <= 1) pos <- c("dummy","dummy")
    if(length(neg) <= 1) neg <- c("dummy","dummy")
    
    toKeep <- c(toKeep, k)
    netReg2[[k]] <- vector("list", 2)
    names(netReg2[[k]]) <- c("pos", "neg")
    netReg2[[k]]$pos <- pos
    netReg2[[k]]$neg <- neg
  }
  netReg2 <- netReg2[toKeep]
  print(length(netReg2))
  return(netReg2)
}

#==========================================================================================
get_regulons <- function(net,corr_matrix,minsize=20){
  regulons=list()
  tnames=c()
  
  for(tfi in rownames(net)){
    tgs <- which(net[tfi,]!=0)
    pos<- which(corr_matrix[tfi,tgs]>0)
    neg<- which(corr_matrix[tfi,tgs]<0)
    r=colnames(net[tfi,net[tfi,]!=0,drop=F])
    if(length(r)>minsize & length(pos)>1 & length(neg)>1) {
      tnames=c(tnames,tfi)
      regulons[[length(regulons)+1]] = r
    }
  }
  names(regulons)=tnames
  return(regulons)
}

#==========================================================================================
#Perform wilcox test on two expression matrices for two cases using identified tfs
perform_wilcox_test <- function(A,B)
{
  genes <- rownames(A);
  wilcox_test_info <- NULL;
  for (i in 1:length(genes))
  {
    temp <- wilcox.test(A[i,],B[i,],exact=FALSE)
    statistic <- temp$statistic
    p_value <- temp$p.value
    if (is.nan(p_value)) { p_value = 1};
    m1 <- mean(A[i,]);
    m2 <- mean(B[i,]);
    FC_mean <- m1-m2;
    wilcox_test_vector <- c(statistic,p_value,m1,m2,FC_mean);
    wilcox_test_info <- rbind(wilcox_test_info,wilcox_test_vector);
  }
  wilcox_test_info <- as.data.frame(wilcox_test_info);
  colnames(wilcox_test_info) <- c("Statistic","Pval","Mean1","Mean2","FC_Mean")
  wilcox_test_info$Statistic <- as.numeric(as.vector(wilcox_test_info$Statistic))
  wilcox_test_info$Pval <- as.numeric(as.vector(wilcox_test_info$Pval))
  wilcox_test_info$Mean1 <- as.numeric(as.vector(wilcox_test_info$Mean1))
  wilcox_test_info$Mean2 <- as.numeric(as.vector(wilcox_test_info$Mean2))
  wilcox_test_info$FC_Mean <- as.numeric(as.vector(wilcox_test_info$FC_Mean))
  wilcox_test_info$Pval <- p.adjust(wilcox_test_info$Pval,method="fdr")
  rownames(wilcox_test_info) <- genes;
  return(wilcox_test_info)
}

#==========================================================================================
#Get the activity matrix
activity_mc <- function(mexp,cormat,tflist=NULL,tau=0.6) {

  mexp.s = mexp
  for(i in 1:nrow(mexp.s)){
    #mexp.s[i,] = (mexp.s[i,] - mean(mexp[i,]))/sd(mexp.s[i,])
    mexp.s[i,] = (mexp.s[i,]-mean(mexp[i,]));
  }
  if (is.null(tflist)) {
    tflist=rownames(cormat)
  }
  actmat = mexp[tflist,]
  actmat[1:length(actmat)]=0
  pb = txtProgressBar(min=1,max=length(tflist),style=3)
  i=1
  for(tfi in tflist) {
    postrg = names(cormat[tfi,cormat[tfi,] > tau])
    negtrg= names(cormat[tfi,cormat[tfi,] < -tau])
    if (length(postrg)>1 & length(negtrg)>1) {
      apos = apply(mexp.s[postrg,,drop=F],2,sum)/length(postrg)
      aneg = apply(mexp.s[negtrg,,drop=F],2,sum)/length(negtrg)
      actmat[tfi,] =  apos - aneg
    } 
    else if (length(postrg)>1 & length(negtrg)<=1)
    {
      actmat[tfi,] = 1;
    }
    else if (length(postrg)<=1 & length(negtrg)>1)
    {
      actmat[tfi,] = -1;
    }
    else
    {
      actmat[tfi,] = 0;
    }
    
    setTxtProgressBar(pb, i)
    i=i+1

  }
  return(actmat)
}

#===========================================================================================
createRegulons <- function(net,ccor,minsize=20)
{
  regulon <- vector("list",nrow(net))
  for(i in 1:nrow(net)){
    tgs <- which(net[i,]!=0)
    pos <- which(ccor[i,tgs]>0)
    neg <- which(ccor[i,tgs]<0)
    if((length(pos)>minsize | length(neg)> minsize) & (length(pos)>1 & length(neg)>1)){
      regulon[[i]] = list(pos=names(pos),neg=names(neg))
    }
    else{
      regulon[[i]] = NA
    }
  }
  names(regulon)<- rownames(net)
  regulon <- regulon[!is.na(regulon)]
  return(regulon)
}

#============================================================================================
massiveGST <- function(rrnk, GGO, alternative = "greater", keepDetails = FALSE, 
                       writeXLS = FALSE, fName = "massiveGST.xls") {
  size <- length(GGO)
  GO <-  intersect(GGO, names(rrnk))
  actualSize <- length(GO)
  result <- data.frame(collection=c("Set"), size, actualSize)
  rnk <- rank(rrnk)
  sumOfRanks <-  sum(rnk[GO])
  n <- length(rnk)
  nx <- actualSize #ny
  ny <- n - nx #nx
  U_stat <- nx * ny + ny * (ny + 1)/2 + sumOfRanks - n * (n + 1)/2
  pod <- U_stat/nx/ny 
  odd <- pod/(1-pod)
  log2_odd <- log2(odd)
  zValue <- U_stat - nx * ny/2
  sigma <- sqrt(nx * ny * (nx + ny + 1)/12)
  correction <- switch(alternative, two.sided = sign(zValue) * 0.5, greater = 0.5, less = -0.5)
  zValue <- (zValue - correction)/sigma
  pValue <- switch(alternative, less = 1 - pnorm(zValue), greater = pnorm(-zValue), two.sided = 2 * pnorm(-abs(zValue)))
  qValue <- p.adjust(pValue, method = "BH") # fdr
  result <- data.frame(result, sumOfRanks, U_stat, pod, odd, log2_odd, zValue, pValue, qValue)
  rowsToRemove <- which(result[, "actualSize"] == 0)
  if(length(rowsToRemove) > 0) result <- result[-rowsToRemove,]
  if(alternative == "two.sided") {
    abs_log2_odd <- abs(result[, "log2_odd"])
    result <- cbind(result, abs_log2_odd)
    positive <- which(result[, "log2_odd"] > 0)
    order_p <- (rank(result[positive, "log2_odd"]) + rank(result[positive, "actualSize"]) + rank(-log10(result[positive, "pValue"])))/3
    names(order_p) <- rownames(result)[positive]
    negative <- which(result[, "log2_odd"] <= 0)
    order_n <- (rank(-result[negative, "log2_odd"]) + rank(result[negative, "actualSize"]) + rank(- log10(result[negative, "pValue"])))/3
    names(order_n) <- rownames(result)[negative]
    ordering <- c(order_p, -order_n)[rownames(result)]
    result <- cbind(result[, 1:12], ordering)
  } else {
    if(alternative == "greater") {
      ordering <- (rank(result[, "log2_odd"]) + rank(result[, "actualSize"]) + rank(-log10(result[, "pValue"])))/3
      names(ordering) <- rownames(result)[rownames(result)]
      result <- cbind(result, ordering)
    } else {
      ordering <- (rank(-result[, "log2_odd"]) + rank(result[, "actualSize"]) + rank(-log10(result[, "pValue"])))/3
      names(ordering) <- rownames(result)[rownames(result)]
      result <- cbind(result, ordering)
    }
  }
  
  ordering <- order(result[, "ordering"], decreasing = TRUE)
  result <- result[ordering,]
  ##########
  if(!keepDetails) {
    colsToRemove <- c("sumOfRanks", "U_stat", "zValue")
    colsToRemove <- which(colnames(result) %in% colsToRemove)
    result <- result[, -colsToRemove]
  }
  
  if(writeXLS) {
    gstTable <- as.data.frame(result)
    require(WriteXLS)
    WriteXLS("gstTable", ExcelFileName = fName, row.names = TRUE)
  }
  invisible(result)
}

massiveExtGst <- function (rankedList, geneSetUp, geneSetDown, minLenGeneSet = 15) 
{
  
  doubleRankedList <- c(rankedList, -rankedList)
  flag <- c(rep(TRUE, length(rankedList)), rep(FALSE, length(rankedList)))
  oorder <- order(doubleRankedList, decreasing = TRUE)
  doubleRankedList <- doubleRankedList[oorder]
  flag <- flag[oorder]
  hits <- (names(doubleRankedList) %in% geneSetUp) & flag
  hits <- hits | ((names(doubleRankedList) %in% geneSetDown) & 
                    !flag)
  names(doubleRankedList)[which(!hits)] <- paste0(names(doubleRankedList)[which(!hits)], 
                                                  1:sum(!hits))
  geneSet <- c(geneSetUp, geneSetDown)
  ans <- massiveGST(doubleRankedList, geneSet,alternative = "two.sided")
  invisible(ans)
}


#=============================================================================================
computeActivity <- function(net,E,regulons,ncores=16){
  M<- apply(E,1,mean)
  E <- E-M
  tf<- intersect(rownames(net),rownames(E))
  targets  <- intersect(colnames(net),rownames(E))
  net <- net[tf,targets]
  E <- E[targets,]
  
  aMat <- matrix(0,nrow=length(regulons),ncol=ncol(E))
  
  for(i in 1:ncol(E)){
    rankedList <- sort(E[,i], decreasing=T)
    #lans <- mclapply(regulon,function(x) mwwExtGST(rankedList = rankedList,geneSetUp = x$pos,
    #                                             geneSetDown = x$neg),mc.cores = ncores)
    lans <- mclapply(regulons,function(x) massiveExtGst(rankedList = rankedList,geneSetUp = x$pos,
                                                       geneSetDown = x$neg),mc.cores = ncores)
    
    nes <- unlist(lapply(lans,function(x) x$pod ))
    activity <- log2(nes/(1-nes))
    aMat[,i] <- activity
  }
  return(aMat)
}


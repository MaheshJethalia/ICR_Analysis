activity = function(mexp,posReg,negReg,tau=0.6) {
  mexp.s = mexp
  #mexp.s = scale(mexp)
  for(i in 1:nrow(mexp.s)){
    mexp.s[i,] = (mexp.s[i,] - mean(mexp[i,]))/sd(mexp.s[i,])
  }
  # remove na
  nass = apply(mexp.s,1,function(x) sum(is.na(x)))
  mexp.s = mexp.s[nass==0,]
  tflist = union(names(posReg),names(negReg))
  #tflist = tflist[tflist %in% rownames(mexp.s)]
  if (length(tflist)==0) {
    stop("no TF names found in datasets")
  }
  actmat = matrix(0,nrow=length(tflist),ncol=ncol(mexp.s))
  rownames(actmat) = tflist
  colnames(actmat) = colnames(mexp)
  pb = txtProgressBar(min=1,max=length(tflist),style=3)
  i=1
  for(tfi in tflist) {
    postrg = posReg[[tfi]]
    postrg = postrg[postrg %in% rownames(mexp.s)]
    
    apos = 1
    if (length(postrg)>0) {
      apos = apply(mexp.s[postrg,,drop=F],2,sum)/length(postrg)
    }
    negtrg = negReg[[tfi]]
    negtrg = negtrg[negtrg %in% rownames(mexp.s)]
    aneg = 1
    if (length(negtrg)>0) {
      aneg = apply(mexp.s[negtrg,,drop=F],2,sum)/length(negtrg)
    }
    actmat[tfi,] =  apos - aneg
    setTxtProgressBar(pb, i)
    i=i+1
    
  }
  return(actmat)
}



ssMRA = function(mexp, regulons,nperm=10000,ncore="all"){
  require(fgsea)
  require(doParallel)
  if(ncore=="all"){
    ncore = parallel:::detectCores()
    doParallel:::registerDoParallel(cores=ncore)
  } else{
    doParallel:::registerDoParallel(cores=ncore)
  }

  mexp.s = mexp
  #mexp.s = scale(mexp)
  for(i in 1:nrow(mexp.s)){
    mexp.s[i,] = (mexp.s[i,] - mean(mexp[i,]))/sd(mexp.s[i,])
  }
  
  m<- max(unlist(lapply(regulons,length)))
  act<-foreach (i=1:ncol(mexp.s), .combine=cbind) %dopar% 
  {
    s<- mexp.s[,i]
    pheno<-sort(s,decreasing=T)
    fgseaRes = fgsea(regulons, pheno, minSize=15, 
                     maxSize=m, nperm=nperm)
    list(fgseaRes) 
    
  }
  ESm<- matrix(NA,nrow=length(regulons),ncol=ncol(mexp.s))
  rownames(ESm)<- names(regulons)
  for(i in 1:length(act)){    
    ESm[act[[i]]$pathway ,i] <- act[[i]]$ES
  }
  NESm<- matrix(NA,nrow=length(regulons),ncol=ncol(mexp.s))
  rownames(NESm)<- names(regulons)
  for(i in 1:length(act)){    
    NESm[act[[i]]$pathway ,i] <- act[[i]]$NES
  }
  PVm<- matrix(NA,nrow=length(regulons),ncol=ncol(mexp.s))
  rownames(PVm)<- names(regulons)
  for(i in 1:length(act)){    
    PVm[act[[i]]$pathway ,i] <- act[[i]]$pval
  }
  colnames(ESm)<-colnames(mexp) 
  colnames(NESm)<-colnames(mexp) 
  colnames(PVm)<-colnames(mexp) 
  ret<-list(ESm,NESm,PVm) 
  names(ret) <- c("ES","NES","PV") 
  return(ret)
}

loadPangliomaExp <- function() {
  e <- load(url("https://www.dropbox.com/s/d3w4eksea16o1e9/gep_1250_agi_rnaseq-dup.RData?raw=1"))
  return(eval(parse(text=e)))
}

loadCasesTable <- function() {
  e <- load(url("https://www.dropbox.com/s/nouain8bh0vfy04/casesTable.RData?raw=1"))
  return(eval(parse(text=e)))
}

loadMechanisticNet <- function() {
  e <- load(url("https://www.dropbox.com/s/c6vag7xyasxup89/me.net.Rdata?raw=1"))
  return(eval(parse(text=e)))
}

loadAracnePangliomaMecNet <- function() {
  e <- load(url("https://www.dropbox.com/s/c2ezk3ehg2kc6nv/mi-aracne-meccnet-tcgaglioma.Rdata?raw=1"))
  return(eval(parse(text=e)))
}


loadPangliomaMC3MAF <- function(){
  e <- load(url("https://www.dropbox.com/s/g95y3r7a26gr5z1/maf_mc3_glio.RData?raw=1"))
  return(eval(parse(text=e)))
}





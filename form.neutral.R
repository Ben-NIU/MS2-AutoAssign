form.neutral<-function(md){
  sp<-strsplit(md, " ")[[1]]
  sap<-sapply(sp, strsplit, split="")
  ## add1 functions to add a number "1" to each 
  add1<-function(x){
    if(is.na(suppressWarnings(as.numeric(x[length(x)])))){
      x[length(x)+1]<-"1"}
    return(x)}
  lap<-lapply(sap, add1)
  element<-lapply(lap, function(x) {paste(x[x %in% c(LETTERS, letters)], collapse = "")})
  num<-lapply(lap, function(x) {as.numeric(paste(x[!x %in% c(LETTERS, letters)], collapse=""))})
  names(num)<-unlist(element)
  g<-as.list(rep(0, 10))
  names(g)<-c("C","H","N","O","S","P","Br","Cl","Si","F")
  gg<-c(g[!names(g) %in% names(num)], num)
  ggg<-list("C"=gg$C, "H"=gg$H,"N"=gg$N,"O"=gg$O,"S"=gg$S,"P"=gg$P,"Br"=gg$Br,"Cl"=gg$Cl,"Si"=gg$Si,"F"=gg$F)
  return(ggg)
}
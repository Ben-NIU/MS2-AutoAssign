Res.to.num<-function(x){
  as.numeric(substr(x, 2, length(strsplit(x, "")[[1]])))
}
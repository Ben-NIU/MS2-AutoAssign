add.Aions<-function(x){
  add.from<-subset(x, Charge==1 & Type=="b" & pos<=2)
  add.from$ms2type<-c("[a1]1+", "[a2]2+")
  add.from$ms2mz<-round(add.from$ms2mz-MonoisotopicMass(formula=list(C=1, O=1)),4)
  add.from$Type<-"a"
  return(rbind(x, add.from))
}
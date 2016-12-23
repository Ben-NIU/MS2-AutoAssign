## This function takes the peptide sequence, modification position and formula and name,and IAA or not to determine the "insilico" y and b series (or c and z ion series) product-ions.
source("FragmentPeptide2.R")
source("form.neutral.R") ## the form.neutral.R takes just the formula of neutral loss, and converts it to the list.
inSilico<-function(pp.seq="", mod.pos=1, mod.fm="", mod.name="X",IAA=TRUE, neutral.fm=""){
  ## this form() function takes the peptide sequence, the position of modification (in number, e.g., 3), and the formula of modification, to generate the final formula of that residue plus the modification.  
  form<-function(mod.fm, pp.seq, mod.pos){
    res<-strsplit(pp.seq, split="")[[1]][mod.pos]
    res.fm<-ConvertPeptide(res, IAA=FALSE)
    if(mod.fm==""){
      gg<-list("C"=res.fm$C, "H"=res.fm$H-2, "N"=res.fm$N, "O"=res.fm$O-1, "S"=res.fm$S, "P"=0,"Br"=0, "Cl"=0, "Si"=0, "F"=0)
    } else {
      sp<-strsplit(mod.fm, " ")[[1]]
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
      gg<-list("C"=gg$C+res.fm$C, "H"=gg$H+res.fm$H-2,"N"=gg$N+res.fm$N,"O"=gg$O+res.fm$O-1,"S"=gg$S+res.fm$S,"P"=gg$P,"Br"=gg$Br,"Cl"=gg$Cl,"Si"=gg$Si,"F"=gg$F)}
    return(gg)
  }
  ## ========== End of form() function.
 if(is.null(neutral.fm)){neutral.fm<-""
 } else {neutral.fm<-neutral.fm}
  
  if(length(strsplit(pp.seq,"")[[1]]) ==1){return(NULL)}
  
  F1<-form(mod.fm, pp.seq, mod.pos)
  F2<-form.neutral(neutral.fm)
  mass<-MonoisotopicMass(formula=F1)   ## this is the mass of that modified residue plus the modification.
  mass.neutral<-MonoisotopicMass(formula=list(C=F1$C-F2$C, H=F1$H-F2$H, N=F1$N-F2$N, O=F1$O-F2$O, S=F1$S-F2$S, P=F1$P-F2$P, Br=F1$Br-F2$Br, Cl=F1$Cl-F2$Cl, Si=F1$Si-F2$Si, F=F1$F-F2$F)) ## this is considering the neutral loss on the same modified residue, refering to the mass of that residue plus the neutral loss.
  
  Len<-strsplit(pp.seq, split="")[[1]]
  Leng<-length(Len)
  if(Leng>mod.pos & mod.pos!=1){
    Mod.pp<-paste(substr(pp.seq, 1, mod.pos-1), "x",substr(pp.seq, mod.pos+1, Leng), sep="")}
  if(Leng>mod.pos & mod.pos==1){
    Mod.pp<-paste("x", substr(pp.seq, 2, Leng), sep="")}
  if(Leng==mod.pos){
    Mod.pp<-paste(substr(pp.seq, 1, Leng-1),"x", sep="")}

  FP2.mod<-FragmentPeptide2(Mod.pp, IAA=IAA, custom=list(code="x", mass=mass))
  FP2.neutral<-FragmentPeptide2(Mod.pp, IAA=IAA, custom = list(code="x", mass=mass.neutral))
  
  FP2.mod<-subset(FP2.mod, select=c("ms2type", "ms2mz"))
  FP2.neutral<-subset(FP2.neutral, select=c("ms2type", "ms2mz"))
  neutral<-FP2.neutral[!FP2.neutral$ms2mz %in% FP2.mod$ms2mz,]
  
  if(nrow(neutral)!=0){
  neutral$ms2type<-paste(neutral$ms2type, "*", sep="")
  FP2.mod<-rbind(FP2.mod, neutral) ## Now combined both dataset with and without the neutral loss, the neutral loss was labeled with "*".
  } else {FP2.mod<-FP2.mod}
  
  charge<-function(x){
    y<-strsplit(as.character(x), "+", fixed = TRUE)[[1]][1]
    yy<-length(strsplit(y, "")[[1]])
    as.numeric(substr(y, yy,yy))
     }
  FP2.mod$Charge<-do.call("c",lapply(FP2.mod$ms2type, charge))
  
  type<-function(x){
    substr(as.character(x), 2,2)
  }
  FP2.mod$Type<-do.call("c", lapply(FP2.mod$ms2type, type))

  pos<-function(x){
   b<-strsplit(strsplit(as.character(x), "]", fixed = TRUE)[[1]][1], "")[[1]]
   as.numeric(paste(b[!is.na(suppressWarnings(as.numeric(b)))],collapse = ""))
  }
  FP2.mod$pos<-do.call("c",lapply(FP2.mod$ms2type, pos))
  
label<-function(x, pp.seq, mod.pos, mod.name){
  if(x[4]=="b" & as.numeric(x[5])>=mod.pos){LB<-paste("+", mod.name, sep="")}
  else if(x[4]=="y" & as.numeric(x[5])>= length(strsplit(pp.seq,"")[[1]]) - mod.pos + 1) {LB<-paste("+", mod.name, sep="")}
  else{LB<-""}
  return(LB)
}
if(mod.fm==""){FP2.mod$label<-""
} else {
  LB<-apply(as.matrix(FP2.mod), 1, label,pp.seq, mod.pos, mod.name)
  FP2.mod$label<-LB}
  return(FP2.mod)
}

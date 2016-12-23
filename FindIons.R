## This is for the Shiny Apps.
## This function takes the scan number of MS2 spectrum, and the peptide sequence, and MS2 tolerance, to determine the y and b series (or c and z ion series) product-ions.
source("FragmentPeptide2.R")
source("form.neutral.R")

FindIons<-function(allPeaks,scan=1234,pp.seq="", mod.pos=1, mod.fm="", mod.name="X", tolerance=50, cutoff=1, IAA=TRUE, neutral.fm=""){

if(is.null(neutral.fm)){neutral.fm<-""
} else {neutral.fm<-neutral.fm}
  

if(is.na(scan) | pp.seq=="") {tb<-NULL
} else {
  if(is.null(tolerance)) {tolerance<-50
  } else {tolerance<-tolerance}
  if(is.null(cutoff) ) {cutoff<-2
  } else {cutoff<-cutoff}
  
  ori<-round(as.data.frame(allPeaks[[scan]]),4)
  names(ori)<-c("mz","int")
  ori$percent<-ori$int*100/max(ori$int)
  ori<-subset(ori, percent>0)
  t<-subset(ori, percent>cutoff)
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
  
  F1<-form(mod.fm, pp.seq, mod.pos)     ## this is the formula of that modified residue plus the modification.
  F2<-form.neutral(neutral.fm)  ## this is just the formula of neutral loss species.
  mass<-MonoisotopicMass(formula=F1)   ## this is the mass of that modified residue plus the modification.
  mass.neutral<-MonoisotopicMass(formula=list(C=F1$C-F2$C, H=F1$H-F2$H, N=F1$N-F2$N, O=F1$O-F2$O, S=F1$S-F2$S, P=F1$P-F2$P, Br=F1$Br-F2$Br, Cl=F1$Cl-F2$Cl, Si=F1$Si-F2$Si, F=F1$F-F2$F)) ## this is considering the neutral loss on the same modified residue, refering to the mass of that residue plus the modification plus the neutral loss.
  
  Len<-strsplit(pp.seq, split="")[[1]]
  Leng<-length(Len)
  if(Leng>mod.pos & mod.pos!=1){
    Mod.pp<-paste(substr(pp.seq, 1, mod.pos-1), "x",substr(pp.seq, mod.pos+1, Leng), sep="")}
  if(Leng>mod.pos & mod.pos==1){
    Mod.pp<-paste("x", substr(pp.seq, 2, Leng), sep="")}
  if(Leng==mod.pos){
    Mod.pp<-paste(substr(pp.seq, 1, Leng-1),"x", sep="")}
  
  FP2<-FragmentPeptide2(pp.seq,IAA=IAA)  ## FP2 refers to the unmodified b/y ion species.
  FP2<-subset(FP2, select=c("ms2type","ms2mz"))
  
  FP2.mod<-FragmentPeptide2(Mod.pp, IAA=IAA, custom=list(code="x", mass=mass)) ## FP2.mod is the b/y ion species when considering the modification on the specific site.
  FP2.mod<-subset(FP2.mod, select=c("ms2type", "ms2mz"))
  
  FP2.neutral<-FragmentPeptide2(Mod.pp, IAA=IAA, custom = list(code="x", mass=mass.neutral)) ##FP2.neutral is the b/y ion series when considering the neutral loss after the modification.
  FP2.neutral<-subset(FP2.neutral, select=c("ms2type", "ms2mz"))
  neutral<-FP2.neutral[!FP2.neutral$ms2mz %in% FP2.mod$ms2mz,]  ## in between FP2.mod and FP2.neutral, only takes the part that shows different, i.e., if neutral loss is "", then these two are the same.
 
  if(nrow(neutral)!=0){
    neutral$ms2type<-paste(neutral$ms2type, "*", sep="")
    FP2.mod<-rbind(FP2.mod, neutral) ## Now combined both dataset with and without the neutral loss, the neutral loss was labeled with "*".
  } else {FP2.mod<-FP2.mod}
  
  ## For the FP2.mod, only retain the m/z that are different from the ones from FP2, so those are ions with modifications.
  # ==================================================================== the following if else identify b/y ions in between experimental spectrum "t" and in silico list FP2/FP2.mod, the if part is when there is no modification (FP2=FP2.mod), the else part is when there is modification (FP2!=FP.
  if(identical(FP2.mod, FP2)){
    ## if this is TRUE, means there are no modification, this is unmodified spectrum.
    tb<-NULL
    for(j in t$mz){
      if( sum(abs(j-FP2$ms2mz)*1E6/j < tolerance) >0 ) {
        id1<-which(abs(j-FP2$ms2mz)*1E6/j < tolerance)
        tb0.unmod<-FP2[id1,]
        tb0<-data.frame("ms2type"=as.character(FP2[id1,]$ms2type), "ms2mz"=FP2[id1,]$ms2mz, "mz"=j, "int"=t$int[which(t$mz==j)])
      } else {
        tb0<-data.frame("ms2type"="","ms2mz"=NA, "mz"=j, "int"=t$int[which(t$mz==j)] )}
      tb<-rbind(tb, tb0)
    }
  } else { ## means there is a modification, the FP2.mod is different from FP2 now.
    FP2.mod<-FP2.mod[!(FP2.mod$ms2mz %in% FP2$ms2mz),]  ## only get the part of FP2.mod that is different from teh FP2.
    FP2.mod$ms2type<-paste(FP2.mod$ms2type, "|", mod.name, sep="")  ##also change the label of those modified ions.
    
    ## find ions (starting from the first peak in ms2 spectrum) that match with the m/z on the in silico list.
    tb<-NULL
    for(j in t$mz){
      if( sum(abs(j-FP2$ms2mz)*1E6/j < tolerance) >0 | sum(abs(j-FP2.mod$ms2mz)*1E6/j < tolerance) >0) {
        id1<-which(abs(j-FP2$ms2mz)*1E6/j < tolerance)
        id2<-which(abs(j-FP2.mod$ms2mz)*1E6/j < tolerance)
        tb0.unmod<-FP2[id1,]
        tb0.mod<-FP2.mod[id2,]
        tb0<-data.frame("ms2type"=c(as.character(FP2[id1,]$ms2type), as.character(FP2.mod[id2,]$ms2type)), "ms2mz"=c(FP2[id1,]$ms2mz, FP2.mod[id2,]$ms2mz), "mz"=j, "int"=t$int[which(t$mz==j)])
      } else {
        tb0<-data.frame("ms2type"="","ms2mz"=NA, "mz"=j, "int"=t$int[which(t$mz==j)] )}
      tb<-rbind(tb, tb0)
    }
  }
  # ========================================================================= the end of this if else, the "tb" contains the identified (matching) m/z.
  
  ## ====== the function BY() tells apart the nature of this fragment, being "none", "b", or "y".
  BY<-function(x) {
    if(as.character(x)==""){
      type<-""
    } else {
      type<-strsplit(as.character(x), split="")[[1]][2]}
    return(type)
  }
  ## ====== End of this BY() function.
  Type<-do.call("c", lapply(tb$ms2type, BY))
  tb$Type<-Type
  ## ===== Creat a function that considers the isotopic peaks of those identified b/y ions, given the charge state of each ion.
  ## === first, create the function that identifies the charge state of each b/y ion.
  Charge<-function(x){
    if(as.character(x)==""){
      cs<-""
    } else {
      part1<-strsplit(as.character(x), split="+", fixed = TRUE)[[1]][1]
      cs<-strsplit(part1, split="")[[1]][length(strsplit(part1, split="")[[1]])]}
    return(cs)}
  charge<-do.call("c", lapply(tb$ms2type, Charge)) 
  tb$Charge<-charge
  ## ======= end of Charge() function.
  ## === second, create the function iso() that takes the original ms2 spectrum (argument 2), and assign the isotopic peaks (A+1, A+2) to those b/y ions.
  
  ## === start of this iso(), the output of this function is the isotopically-labeled b/y ions dataframe.
  iso<-function(ms2){## here the ms2 is the "tb"
    byion<-subset(ms2, Charge!="")$mz
    for(u in byion){
      e<-which(ms2$mz==u)
      if(e+2<=nrow(ms2)){  
        if(abs((ms2$mz[e+1]-ms2$mz[e])*as.numeric(ms2$Charge[e])-1) < 0.01  & ms2$Type[e+1]==""){
          ms2$Type[e+1]<-ms2$Type[e]
          if(abs((ms2$mz[e+2]-ms2$mz[e+1])*as.numeric(ms2$Charge[e])-1) < 0.01 & ms2$Type[e+2]==""){
            ms2$Type[e+2]<-ms2$Type[e+1]}}
      }
    }
    return(ms2)}
  ## === end of this iso().
  
  tb<-suppressWarnings(iso(tb))
  ## ==== this is to make sure that tb contains the whole spectrum.
 empty<-data.frame(ms2type="",ms2mz=NA, mz=ori$mz, int=ori$int, Type="",Charge="")
 tb<-rbind(empty[!empty$mz %in% tb$mz,], tb)
  
  }
  return(tb)
}

## This is for the Shiny Apps.
## This function takes the scan number of MS2 spectrum, and the peptide sequence, and MS2 tolerance, to determine the y and b series (or c and z ion series) product-ions.
source("FragmentPeptide2.R")
source("form.neutral.R")
source("inSilico.R")
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

   
    FP2.mod<-inSilico(pp.seq, mod.pos, mod.fm, mod.name, IAA=IAA, neutral.fm)
    FP2.unMod<-inSilico(pp.seq, mod.pos, mod.fm="", mod.name="", IAA=IAA, neutral.fm="")
    
  
    ## For the FP2.mod, only retain the m/z that are different from the ones from FP2, so those are ions with modifications.
    # ==================================================================== the following if else identify b/y ions in between experimental spectrum "t" and in silico list FP2/FP2.mod, the if part is when there is no modification (FP2=FP2.mod), the else part is when there is modification (FP2!=FP.
    if(identical(FP2.mod, FP2.unMod)){
      ## if this is TRUE, means there are no modification, this is unmodified spectrum.
      tb<-NULL
      for(j in t$mz){
        if( sum(abs(j-FP2.unMod$ms2mz)*1E6/j < tolerance) >0 ) {
          id1<-which(abs(j-FP2.unMod$ms2mz)*1E6/j < tolerance)
          tb0.unmod<-FP2.unMod[id1,]
          tb0<-data.frame("ms2type"=as.character(FP2.unMod[id1,]$ms2type), "ms2mz"=FP2.unMod[id1,]$ms2mz, "mz"=j, "int"=t$int[which(t$mz==j)])
        } else {
          tb0<-data.frame("ms2type"="","ms2mz"=NA, "mz"=j, "int"=t$int[which(t$mz==j)] )}
        tb<-rbind(tb, tb0)
      }
    } else { ## means there is a modification, the FP2.mod is different from FP2.unMod now.
      FP2.mod<-FP2.mod[!(FP2.mod$ms2mz %in% FP2.unMod$ms2mz),]  ## only get the part of FP2.mod that is different from teh FP2.unMod.
      FP2.mod$ms2type<-paste(FP2.mod$ms2type, "|", mod.name, sep="")  ##also change the label of those modified ions.
      
      ## find ions (starting from the first peak in ms2 spectrum) that match with the m/z on the in silico list.
      tb<-NULL
      for(j in t$mz){
        if( sum(abs(j-FP2.unMod$ms2mz)*1E6/j < tolerance) >0 | sum(abs(j-FP2.mod$ms2mz)*1E6/j < tolerance) >0) {
          id1<-which(abs(j-FP2.unMod$ms2mz)*1E6/j < tolerance)
          id2<-which(abs(j-FP2.mod$ms2mz)*1E6/j < tolerance)
          tb0.unmod<-FP2.unMod[id1,]
          tb0.mod<-FP2.mod[id2,]
          tb0<-data.frame("ms2type"=c(as.character(FP2.unMod[id1,]$ms2type), as.character(FP2.mod[id2,]$ms2type)), "ms2mz"=c(FP2[id1,]$ms2mz, FP2.mod[id2,]$ms2mz), "mz"=j, "int"=t$int[which(t$mz==j)])
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

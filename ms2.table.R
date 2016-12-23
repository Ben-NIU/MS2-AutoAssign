ms2.table<-function(from.FindIons, annotated){
 if(annotated=="NO"){FI<-NULL
 } else { ## when annotated=="YES".
   if(is.null(from.FindIons)) {FI<-NULL
       } else {
       FI<-subset(from.FindIons, ms2type!="")
        if(nrow(FI)==0){FI<-NULL
                } else {
                 pos<-function(x){
                 b<-strsplit(strsplit(as.character(x), "]", fixed = TRUE)[[1]][1], "")[[1]]
                 as.numeric(paste(b[!is.na(suppressWarnings(as.numeric(b)))],collapse = ""))
                   }
                   FI$pos<-do.call("c",lapply(FI$ms2type, pos))

                  FI<-subset(FI[order(FI$Type, FI$pos, FI$Charge),], select=c("ms2type","mz", "int"))
                  FI$mz<-round(FI$mz, 4)
                  FI$int<-round(FI$int,4)
                       }}
  }
  return(FI)
}
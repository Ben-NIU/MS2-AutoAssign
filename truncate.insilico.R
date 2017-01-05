truncate.insilico<-function(fulllist,charge, type, NT){
  if(is.null(fulllist)){list2<-NULL
  } else {
    
  with.nt<-function(x){
    "*" %in% strsplit(as.character(x),"")[[1]]
  }
  
 fulllist$withNT<-do.call("c", lapply(fulllist$ms2type, with.nt))
  list1<-rbind(subset(fulllist, Charge==charge[1]), subset(fulllist, Charge==charge[2]), subset(fulllist, Charge==charge[3]))
  list2<-rbind(subset(list1, Type==type[1]), subset(list1, Type==type[2]), subset(list1, Type==type[3]))
  list2<-list2[order(list2$Type, list2$pos),]
 if(NT==TRUE){list2<-list2
 } else {list2<-subset(list2, withNT==FALSE)}
  }
  
  return(list2)
}
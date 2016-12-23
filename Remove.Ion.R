Remove.Ion<-function(findIon, tableIon, row.sel){ ## this is to remove the labels (ms2type) of those selected ions.
  if(is.null(findIon) | is.null(tableIon)){
    return(NULL) } else {
  findIon$ms2type[findIon$mz %in% tableIon$mz[row.sel]]<-""
  return(findIon)
    }
}
  
  
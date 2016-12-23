## This is for the Shiny Apps
## ============== Below function Plot.FindIons() will plot the annotated MS2 spectrum.
Plot.FindIons<-function(from.FindIons,allPeaks,scan, pp.seq, left, right, force=4, annotated="NO", alpha,adjust.x, adjust.y, label.size=4, peakwd=0.8){

if(is.null(from.FindIons) & is.na(scan)) {return(NULL)}
  
else if(is.null(from.FindIons) & !is.na(scan)) { ## indicate that pp.seq=="" in this case.
   dt<-as.data.frame(allPeaks[[scan]])
   names(dt)<-c("mz","int")
   ggplot(dt, aes(x=mz, xend=mz, y=0, yend=int)) + geom_segment(size=peakwd)+theme(panel.border=element_rect(fill=NA, size=1, linetype="solid", color="black"), panel.background=element_blank(), panel.grid.major=element_line(color="grey95")) + theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=20), axis.title.x=element_text(face="bold", size=20), axis.title.y=element_text(face="bold", size=20), plot.subtitle = element_text(hjust=0.5, size=15)) + labs(x="m/z",y="Relative Intensity", subtitle=paste("scan #", scan, sep="")) + theme(legend.title=element_text(color="black",size=20),legend.text=element_text(size=20), legend.position="none") +scale_y_continuous(name="Relative Intensity", breaks=max(dt$int)/10*(c(0,5,10)), labels=c(0,50,100), expand=c(0,0)) + coord_cartesian(ylim=c(0, max(dt$int)*1.1)) + xlim(left, right)
} 
else{ ## indicate that both pp.seq and scan exist now.
  dt<-from.FindIons
  if(annotated=="NO"){
    ggplot(dt, aes(x=mz, xend=mz, y=0, yend=int, color=Type)) + geom_segment(size=peakwd)+theme(panel.border=element_rect(fill=NA, size=1, linetype="solid", color="black"), panel.background=element_blank(), panel.grid.major=element_line(color="grey95")) + theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=20), axis.title.x=element_text(face="bold", size=20), axis.title.y=element_text(face="bold", size=20), plot.subtitle = element_text(hjust=0.5, size=15)) + labs(x="m/z",y="Relative Intensity", subtitle=pp.seq) + scale_color_manual(values=c("black","#0072B2","#CC79A7"), breaks=c("","b","y"), labels=c("unidentified","b-ion","y-ion"), name="") + theme(legend.title=element_text(color="black",size=20),legend.text=element_text(size=20), legend.position="none") +scale_y_continuous(name="Relative Intensity", breaks=max(dt$int)/10*(c(0,5,10)), labels=c(0,50,100), expand=c(0,0)) + coord_cartesian(ylim=c(0, max(dt$int)*1.1)) + xlim(left, right)
  } 
   
  else if(annotated=="YES" & nrow(subset(subset(dt,mz>=left & mz<=right), ms2type!=""))==0){
  ggplot(dt, aes(x=mz, xend=mz, y=0, yend=int, color=Type)) + geom_segment(size=peakwd)+theme(panel.border=element_rect(fill=NA, size=1, linetype="solid", color="black"), panel.background=element_blank(), panel.grid.major=element_line(color="grey95")) + theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=20), axis.title.x=element_text(face="bold", size=20), axis.title.y=element_text(face="bold", size=20), plot.subtitle = element_text(hjust=0.5, size=15)) + labs(x="m/z",y="Relative Intensity", subtitle=pp.seq) + scale_color_manual(values=c("black","#0072B2","#CC79A7"), breaks=c("","b","y"), labels=c("unidentified","b-ion","y-ion"), name="") + theme(legend.title=element_text(color="black",size=20),legend.text=element_text(size=20), legend.position="none") +scale_y_continuous(name="Relative Intensity", breaks=max(dt$int)/10*(c(0,5,10)), labels=c(0,50,100), expand=c(0,0)) + coord_cartesian(ylim=c(0, max(dt$int)*1.1)) + xlim(left, right)
  } else {
  ggplot(dt, aes(x=mz, xend=mz, y=0, yend=int, color=Type)) + geom_segment(size=peakwd)+theme(panel.border=element_rect(fill=NA, size=1, linetype="solid", color="black"), panel.background=element_blank(), panel.grid.major=element_line(color="grey95")) + theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=20), axis.title.x=element_text(face="bold", size=20), axis.title.y=element_text(face="bold", size=20), plot.subtitle = element_text(hjust=0.5, size=15)) + labs(x="m/z",y="Relative Intensity", subtitle=pp.seq) + scale_color_manual(values=c("black","#0072B2","#CC79A7"), breaks=c("","b","y"), labels=c("unidentified","b-ion","y-ion"), name="") + theme(legend.title=element_text(color="black",size=20),legend.text=element_text(size=20), legend.position="none") +scale_y_continuous(name="Relative Intensity", breaks=max(dt$int)/10*(c(0,5,10)), labels=c(0,50,100), expand=c(0,0)) + coord_cartesian(ylim=c(0, max(dt$int)*1.1)) + xlim(left, right)+geom_text_repel(aes(y=int,label=ms2type,color=Type), angle=90, show.legend = FALSE, nudge_x = adjust.x, nudge_y = adjust.y, force=force, segment.alpha = alpha, size=label.size)
  }
}
}




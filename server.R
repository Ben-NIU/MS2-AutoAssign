library(OrgMassSpecR)
library(mzR)
library(ggplot2)
library(ggrepel) 
source("inSilico.R")
source("Res.to.num.R")
source("truncate.insilico.R")
source("FindIons.R")
source("Plot.FindIons.R")
source("ms2.table.R")
source("Remove.Ion.R")
## this is to get rid of the upper limit for file uploading.
options(shiny.maxRequestSize=-1)

shinyServer(function(input, output) {
  
  residues<-reactive({
    validate(
      need(input$ppseq !="", "Note: two or more residues required!"))
    paste(strsplit(input$ppseq,"")[[1]], 1:length(strsplit(input$ppseq,"")[[1]]), sep="") })
  
  output$Mod.pos<-renderMenu({
    if(residues()[1]=="") return()
    if(residues()[1]!=""){
    selectInput("Mod.pos", label="site modified", choices = residues(),selected=NULL)}
  })
  
  output$Neutral<-renderUI({
    if(input$NT=="No"){return()}
    if(input$NT=="Yes"){textInput("Neutral",label="formula of neutral loss", value="", placeholder = "C H2 N2")}
  })
  
  num<-reactive({Res.to.num(input$Mod.pos)})
  
  full.list<-reactive({
    validate(
      need(input$ppseq !="", ""))
    inSilico(input$ppseq, num(), input$mod.fml, input$mod.name, input$IAA, input$Neutral)
  })
  
  observeEvent(input$bt,{
  output$tb2<-DT::renderDataTable({
      isolate(truncate.insilico(full.list(), charge=input$charge, type=input$type, NT=input$withNT)[,c(1,2,6)] )}, rownames=FALSE,options=list(
        lengthMenu=list(c(10,25,-1),c("10","25","All")),
        pageLength=25))})
      
### Read the uploaded .mzXML file.
  Peaks<-reactive({
    validate(
      need(input$mzXML$datapath != "", "")
    )
    peaks(openMSfile(input$mzXML[['datapath']])) })
  
  HD<-reactive({
    validate(
      need(input$mzXML$datapath != "", "")
    )
    header(openMSfile(input$mzXML[['datapath']])) })
  
  findIons<-reactive({ FindIons(Peaks(), input$scan, input$ppseq, num(), input$mod.fml, input$mod.name, input$tol, input$cutoff, input$IAA, input$Neutral)})
  
  tableIons<-reactive({ms2.table(findIons(), input$annotated)})
  
## =====following are some dynamic UIs, depending on the status of "annotated".
     output$tol<-renderUI({ ## if annotated is selected as "NO", then no UI will be rendered
        if(input$annotated=="NO") return()
        if(input$annotated=="YES") {numericInput("tol", label="tolerance in ppm", value=50, min=3, max=500, step=1)}
        })
      output$cutoff<-renderUI({
        if(input$annotated=="NO") return()
        if(input$annotated=="YES") {sliderInput("cutoff", label="intensity cutoff %", min=0, max=35, value=5, step=1)}
      })
      output$label.size<-renderUI({
        if(input$annotated=="NO") return()
        if(input$annotated=="YES") {numericInput("label.size", label="labeling size", value=4, min=1, max=10, step=1)}
      })
      output$force<-renderUI({
        if(input$annotated=="NO") return()
        if(input$annotated=="YES") {numericInput("force", label="Repulsion", value=3, min=1, max=15, step=1)}
      })
      output$alpha<-renderUI({
        if(input$annotated=="NO") return()
        if(input$annotated=="YES") {numericInput("alpha", label="Transparency", value=0.8, min=0, max=1, step=0.1)}
      })
      output$adjust.x<-renderUI({
        if(input$annotated=="NO") return()
        if(input$annotated=="YES") {numericInput("adjust.x", label="nudge-x", value=0, step=1)}
      })
      output$adjust.y<-renderUI({
        if(input$annotated=="NO") return()
        if(input$annotated=="YES") {numericInput("adjust.y", label="nudge-y", value=5, step=1)}
      })
      
## ====== End of annotated-dependent dynamic UIs.

    output$remove<-renderUI({
      if(input$`remove?`=="No") return()
      if(input$`remove?`=="Yes") {actionButton("remove", label = "Click to update!", icon=icon("refresh"))}
    })  
  
    output$help1<-renderUI({
      if(input$`remove?`=="No") return()
      if(input$`remove?`=="Yes") {helpText(span("How: Select unwanted ions/annotations from Table above and click this button", style="font-family:'calibri';color:grey50; font-size:11pt"))}
    }) 
## ====== the above two dynamic UI are for the options of selective removal ions.
    
    output$RTstat<-renderValueBox({
      if(is.null(input$scan) | is.na(input$scan)) {
        valueBox( value=tags$p("...", style="font-size: 70%"), "Retention Time (min)", icon=icon("clock-o"), color = "purple")
      } else {
      valueBox( value=tags$p(round(HD()$retentionTime[input$scan]/60,2),style="font-size: 70%"), "Retention Time (min)",icon=icon("clock-o"), color="purple")}
    })
    
    output$PrecursorMS<-renderValueBox({
      if(is.null(input$scan) | is.na(input$scan)) {
        valueBox( value=tags$p("...", style="font-size: 70%"), "Precursor m/z", icon=icon("qq"), color = "olive")
       } else {
      valueBox( value=tags$p(round(HD()$precursorMZ[input$scan], 4), style="font-size: 70%"), "Precursor m/z", icon=icon("qq"), color = "olive")}
    })
    
    output$PrecursorCS<-renderValueBox({
      if(is.null(input$scan) | is.na(input$scan)) {
        valueBox( value=tags$p("...", style="font-size: 70%"), "Charge state", icon=icon("lightbulb-o"), color = "yellow")
      } else {
        valueBox( value=tags$p(paste("+ ",round(HD()$precursorCharge[input$scan], 4), sep=""), style="font-size: 70%"), "Charge state", icon=icon("lightbulb-o"), color = "yellow")}
    })
    
 output$tab.spectrum<-renderPlot({
        Plot.FindIons(findIons(), Peaks(),scan=input$scan, input$ppseq, input$left, input$right, input$force, input$annotated, input$alpha, input$adjust.x, input$adjust.y, input$label.size, input$peakwidth)
      })

 output$tab.table<-DT::renderDataTable({
  tableIons()}, rownames=FALSE,options=list(paging=FALSE))
 
 row.selec<-reactive({as.numeric(input$tab.table_rows_selected)})
       
 removed.findIons<-reactive({Remove.Ion(findIons(), tableIons(), row.selec())})
 
 observeEvent(input$remove,{
   output$new.spectrum<-renderPlot({
     isolate(Plot.FindIons(removed.findIons(), Peaks(),scan=input$scan, input$ppseq,input$left, input$right, input$force, input$annotated, input$alpha, input$adjust.x, input$adjust.y, input$label.size, input$peakwidth)) 
   }) })
 ## ====== this is for plot downloading using ggsave().
 output$down.spec<-downloadHandler(filename = function(){paste("MS2-",Sys.Date(),".png", sep="")}, content = function(file){
   device<-function(...,width=width, height=height) grDevices::png(..., width = width, height = height, res=input$Res, units = "in")
   ggsave(file,plot=isolate(Plot.FindIons(removed.findIons(), Peaks(),scan=input$scan, input$ppseq,input$left, input$right, input$force, input$annotated, input$alpha, input$adjust.x, input$adjust.y, input$label.size, input$peakwidth)),device ="png", dpi = input$Res)}, contentType = "image/png")  
 ## ====== thsi is for table downloading using write.csv().
 output$down.table<-downloadHandler(filename = function(){paste("PeakList-", Sys.Date(), ".csv", sep="")}, content = function(file){write.csv(tableIons(), file)})
  
}
)
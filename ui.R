library(shinydashboard)

dashboardPage(
  dashboardHeader(title="MS/MS Spectrum Analysis & Assignment", titleWidth = 450),
  dashboardSidebar(
    sidebarMenu(
    
      sidebarSearchForm("ppseq", "bt", label="Input peptide sequence...", icon=icon("search")),
      radioButtons("IAA", label="Carbamiodomethyl Cys ?", choices = c("Yes"=TRUE,"No"=FALSE), selected = TRUE, inline = TRUE),
      menuItem("Instructions", tabName ="welcome", icon=icon("info-circle"), selected = TRUE ),
      menuItem("MS/MS", tabName = "ms2",icon=icon("gears"),
               menuItemOutput("Mod.pos"),
               textInput("mod.fml", label="formula of labeling", value="", placeholder = "C2 H3 N O"),
               textInput("mod.name", label="name your labeling", value="", placeholder = "GEE"),
               radioButtons("NT", label="Search Neutral Loss?", choices = c("Yes", "No"), selected = "No", inline = TRUE),
               uiOutput("Neutral"),
               menuSubItem("in silico Fragments",tabName = "insilico", icon=icon("th")),
               menuSubItem("Peaks Matching", tabName = "pkmatch", icon=icon("bar-chart"))),
      menuItem("Precursor MS", icon=icon("gear"), badgeLabel = "website",href="https://benniu-wustl.shinyapps.io/proteomics-toolkit/")
 
  )
  ),
  dashboardBody(

    
    tabItems(
      
      tabItem(tabName = "welcome",
              fluidRow(
              box(title="About this website", status="primary", solidHeader = TRUE, collapsible = FALSE, width = 12,
                  p(span(strong("General description"), style="font-family:'calibri';color:black; font-size:12pt")),
                  p(span("This webpage is designed to assist MS/MS computation and analysis, its functionalities include",span(em("in silico"),style="font-family:'times';color:grey50; font-size:11pt"), span("generation of b/y-ions from CID fragmentation, spectra visualization, automatic MS/MS peak assignments, MS/MS neutral loss peak matching, selective annotation, spectrum output (with user-defined resolution), etc.",style="color:grey50; font-size:11pt"),style="color:grey50; font-size:11pt"))
              )),
              fluidRow(
                column(width=7,
                box(title="Brief Tutorial", status = "primary", solidHeader = TRUE, collapsible = TRUE, width=NULL,
                  p(span(strong("MS/MS"),style="font-family:'calibri';color:black; font-size:12pt")),
                  p(span("This tab contains two sub-items, i.e., in silico Fragments and Peak Matching.", style="color:grey50; font-size:11pt")),
                  p(span(strong("in silico Fragments:"), style="color:black; font-size:11pt")),
                  p(span("Input peptide sequence from top of sidebar. For labeled peptide, need specify the chemical formula (e.g., C2 H3 N O, the space in between elements are required!) and name of labeling. Neutral loss can be considered by inputting the net composition of neutral loss.", style="color:grey50; font-size:10.5pt")),
                  p(span(strong("Peaks Matching:"), style="color:black; font-size:11pt")),
                  p(span("This tab focuses on MS/MS spectrum display, annotation, and output. First upload file in .mzXML format.", a("Use ProteoWizard to convert raw data to mxXML.", href="http://proteowizard.sourceforge.net/downloads.shtml", target="_blank"), "The spectrum with given scan number will be shown in MS/MS-panel, one can adjust the displayed m/z range and MS peak width in Spectrum Settings. When there is no peptide sequence input, the scan number is displayed as spectral title; however, a peptide sequence is needed in order to assign the peaks.",style="color:grey50; font-size:10.5pt")),
                  p(span("Within Annotation Settings, the intensity cutoff% determines the shreshold of relative intensities for peaks to be annotated. Tolerance in ppm refers to the tolerance of m/z difference between in silico and experimental results. Repulsion works in a way as if all the annotated labels would repel, so that they won't stack onto each other, the higher this value is, the more repulsion would be.", style="color:grey50; font-size:10.5pt"))
                  
                  )
              ),
        
                column(width=5,
                  box(title="About Us", status = "primary", width=NULL,   
                      hr(),
                      img(src="wustl-logo.png", width=90, height=90, align="right"), 
                      br(),
                       p(span(a("Washington University in St. Louis", href="https://wustl.edu", target="_blank"), style="font-family:'times';color:black; font-size:10.5pt")),
                       p(span(a("Gross Lab: NIH/NCRR Mass Spectrometry Resource", href="http://msr.dom.wustl.edu/michael-l-gross-ph-d/", target="_blank"), style="font-family:'times';color:black; font-size:10.5pt"))
                      
                      
                       ),
            
                 fluidRow(
                   column(width=8,
                  box(title="About the Author", width=NULL, status="warning",
                      hr(),
                      p(span(em("Mr. Ben NIU"), style="font-family:'chalkboard';color:black; font-size:10.5pt")),
                      p(span(em("email: benniu720@gmail.com"),style="font-family:'chalkboard';color:black; font-size:10.5pt")),
                      p(span("Ph.D. Candidate, Gross Lab", style="font-family:'times';color:black; font-size:10.5pt")),
                      p(span("Washington University in St. Louis", style="font-family:'times';color:black; font-size:10.5pt"))
                      
                      )))
                ))
              ),
      
      
      
      
      
      
      
      
      
      tabItem(tabName = "insilico",
                fluidRow(
                  box(title="in silico fragments", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                      DT::dataTableOutput("tb2")),
                  box(title="Display Settings", status="primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                      selectInput("charge", label="Charge state to display", choices = c("+1"="1","+2"="2","+3"="3"), selected = c("1","2"), multiple = TRUE),
                      fluidRow(
                        column(5,checkboxGroupInput("type", label="Ion series to display", choices = c("y-ion"="y","b-ion"="b", "a-ion"="a"), selected =c("y","b"), inline = TRUE  )),
                        column(5, checkboxInput("withNT", label="Show neutral loss?", value=FALSE))
                      ))
              
              )),
      tabItem(tabName = "pkmatch",
          fluidRow(
            column(width=4,
            box(title="Upload the .mzXML file here", width=NULL,background = "light-blue", collapsible = TRUE,
                fileInput("mzXML", label=NULL),
                hr(),
                numericInput("scan", label="Input spectrum id.",value="" )),
            box(title="Spectrum Settings", collapsible = TRUE, width=NULL, solidHeader = TRUE, status = "primary",
                fluidRow(
                  column(width=6, numericInput("left", label="lower m/z", value=80)),
                  column(width=6, numericInput("right", label="upper m/z", value=2200))
                ),
                fluidRow(
                  column(width=12, sliderInput("peakwidth", label="Line width of spectrum", min=0.3, max=1, value=0.8, step=0.1))
                )
            )
            ),
            column(width=8,
            box(title="Annotation Settings", collapsible = TRUE, width = NULL,solidHeader = TRUE, status = "primary", 
                fluidRow(
                  column(5,radioButtons("annotated", label="Annotated spectrum?", choices = c("NO","YES"), selected = "NO", inline = TRUE))
                ),
                hr(),
                fluidRow(
                  column(width=3, uiOutput("tol")),
                  column(width=6, uiOutput("cutoff")),
                  column(width=3, uiOutput("label.size"))
                ),
                fluidRow(
                  column(width=3, uiOutput("force")),
                  column(width=3, uiOutput("alpha")),
                  column(width=3, uiOutput("adjust.x")),
                  column(width=3, uiOutput("adjust.y"))
                )),
            box(title="Spectral Stats", collapsible = TRUE, width=NULL,status="primary",
            fluidRow(
           valueBoxOutput("RTstat"),
           valueBoxOutput("PrecursorMS"),
           valueBoxOutput("PrecursorCS")
            ))
      
          )),
          fluidRow(
            tabBox(width=12,side="right", id="tabbox", selected = "Spectrum", title = span("MS/MS", style="font-family:'calibri';color:black; font-size:12pt"),
              tabPanel("Table", icon = icon("table"),
              DT::dataTableOutput("tab.table")),
              tabPanel("Spectrum", icon=icon("area-chart"),
              plotOutput("tab.spectrum"))
             
            
            )
         
    
          ),
          fluidRow(
            box(title="Options", width=4, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,status = "primary",
    
                radioButtons("remove?", label="Selectively remove peak annotations?", choices = c("Yes","No"), selected="No", inline=TRUE),
                uiOutput("help1"),
                uiOutput("remove")),
            box(title="Spectrum Output", width=4, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,status="warning",
                helpText(span("Click download button to output the High-Res annotated MS/MS spectrum", style="font-family:'calibri';color:grey50; font-size:11pt")),
               
              fluidRow(
               
                column(5, span("Resolution (dpi)", style="font-family:'calibri';color:black; font-size:11pt"),
                       numericInput("Res", label = NULL, value=300)),
                
                column(6, span("Here", style="font-family:'calibri';color:white; font-size:11pt"),       
                       downloadButton("down.spec", label="Download"))
                )),
            
            box(title="Table Output", width=4, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning",
                helpText(span("Click Output Table to download the table for those identified MS/MS peaks", style="font-family:'calibri';color:grey50; font-size:11pt")),
                downloadButton("down.table", label="Output Table"))
                
                
                ),
          fluidRow(
            box(title="MS/MS with selective removal of annotations", width=12,solidHeader = TRUE, collapsible = TRUE,collapsed = TRUE,status="primary",
               plotOutput("new.spectrum")
                ))
          )
              
              )
    
    
  )
)


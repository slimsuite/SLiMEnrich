#*********************************************************************************************************
#*********************************************************************************************************
# Short Linear Motif Enrichment Analysis App (SLiMEnrich)
# Developer: **Sobia Idrees**
# Version: 1.1.1
# Description: SLiMEnrich predicts Domain Motif Interactions (DMIs) from Protein-Protein Interaction (PPI) data and analyzes enrichment through permutation test.
#*********************************************************************************************************
#*********************************************************************************************************
##############################
#Version History
##############################
#V1.0.1 - Added code for checking whether packages installed. (Removes manual step)
#V1.0.2 - Better naming conventions in code
#V1.0.3 - Added titles/captions to data tables (uploaded files).
#       - Improved summary bar chart (used plotly), 
#       - Improved histogram module (removed separate window option for plot, added width/height input option in settings of histogram to download plot as png file). 
#V1.0.4 - Checks whether any of the random files is missing and creates if not present.
#V1.0.5 - Added a new tab to show distribution of ELMs in the predicted DMI dataset in tabular as well as in interactive view.
#V1.0.7 - File headers to lowercase for consistency
#V1.0.8 - Auto loading example dataset
#V1.0.9 - Reads SLiMProb REST server output through Job Id.
#V1.1.0 - Added new tab to show distribution of Domains in the predicted DMI dataset.
#V1.1.1 - New FDR calculation
#V1.2.0 - Uses known and predicted ELM information. Predicts DMIs based on domains as well as based on proteins.
##############################
#SLiMEnrich program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# SLiMEnrich program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################
#Required Libraries
##############################
package_names = c("shiny", "ggplot2", "colourpicker", "shinyBS", "shinythemes", "DT", "shinyjs", "visNetwork", "igraph","markdown","plotly", "plyr", "shinyWidgets")
for(package_name in package_names) 
{ 
  library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE) 
} 

##############################
#GUI of the App
##############################
#navbar page with sidebar layout along with tabsets
ui <- shinyUI(navbarPage(div(id= "title", ("SLiMEnrich")),windowTitle = "SLiMEnrich", tabPanel("Domain-Motif Interactions", tags$head(
  tags$style(HTML("
                  .shiny-output-error-validation {
                  color: red;
                  font-size: 18px;
                  font-style: italic;
                  font-weight: bold;
                  -webkit-animation: mymove 5s infinite; /* Chrome, Safari, Opera */
                  animation: mymove 5s infinite;
                  }
                  @-webkit-keyframes mymove {
                  50% {color: black;}
                  }
                  .shiny-notification {
                  font-weight: bold;
                  bottom: calc(0%)
                  width: calc(100%)
                  }
                  #note{
                  color: red;
                  background-color: white;
                  font-size: 10;
                  font-weight: bold;
                  
                  }
                  #update{
                  color: white;
                  background-color: black;
                  font-weight: bold;
                  font-size: 10;
                  }
                  
                  " ))
  ),
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      fileInput("PPI","Select Interaction file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
      
      
      div(id = "slimrun", textInput("SLiMRun", label = "", value = "")),
      actionButton("run", "Run", width = "100px"),
      hr(),
      
      prettyRadioButtons(inputId = "options",
                         label = "Motif Prediction Type", icon = icon("check"),
                         choices = c("Known SLiMs" = "true", "Predicted SLiMs" = "pred"),
                         animation = "pulse", status = "info"),
      hr(),
      div(class="pretty p-switch p-fill",prettyCheckbox("withPro", label = tags$b("Protein Motif Interaction (PMI)"), value = FALSE, status = "info", 
                                                        
                                                        animation = "pulse")),
      hr(),
      div(id="fileuploads",prettyCheckbox("uploadmotifs",label = tags$b("Upload Files"), value = FALSE, status = "warning",
                                          icon = icon("check"),
                                          animation = "pulse")),
      hr(),
      
      div(id="uploadmotif",  fileInput("domain","Select Domain file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
          fileInput("MotifDomain","Select Motif-Domain file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
          fileInput("Motif","Select SLiM prediction file (e.g. SLiMProb)",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
          prettyCheckbox("SLiMrunid", label = "Provide SLiMProb Job ID", status = "default",
                         icon = icon("check"),
                         animation = "pulse")),
      div (id = "note", "Note: To analyze example dataset, press 'Run' without uploading any files"),
      hr(),
      div (id = "update", "Last updated: 29-Jun-2018")
    ),
    
    # MainPanel
    mainPanel(
      #Creates a seperate window (pop up window)
      bsModal("DisE", "ELM Distribution", "godis", size = "large", plotlyOutput("diselmchart")),
      bsModal("DidsE", "Domain Distribution", "godisd", size = "large", plotlyOutput("disdomchart")),
      #Tab view
      tabsetPanel(type="tabs",
                  tabPanel("Uploaded Data",
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%", "50%", "50%"), DT::dataTableOutput("udata2"), DT::dataTableOutput("udata")), DT::dataTableOutput("udata4"), DT::dataTableOutput("udata3")
                           )
                  ),
                  
                  tabPanel("Potential DMIs",
                           DT::dataTableOutput("data"),
                           tags$hr(),
                           downloadButton('downloadDMI', 'Download')
                  ),
                  
                  tabPanel("Predicted DMIs", DT::dataTableOutput("PredDMIs"),tags$hr(),downloadButton('downloadpredDMI', 'Download')),
                  
                  tabPanel("Statistics", fluidRow(
                    splitLayout(cellWidths = c("75%", "25%"), plotlyOutput("plotbar"))
                  )),
                  tabPanel("Histogram", fluidRow(
                    splitLayout(cellWidths = c("50%", "50%"), plotOutput("histogram"), htmlOutput("summary"))),
                    tags$hr(),
                    div(id="txtbox",actionButton("setting", "Settings")),
                    div(id="txtbox",downloadButton("downloadPlot", "Download")),
                    
                    div(id="settings", sliderInput("bins", 
                                                   "Number of bins",
                                                   min= 1,
                                                   max = 200,
                                                   value = 30),
                        tags$hr(),
                        tags$h4(tags$strong("Select labels")),
                        
                        checkboxInput("barlabel", label="Bar Labels", value = FALSE, width = NULL),
                        div(id="txtbox", textInput("text3", label = "Main title", value = "Distribution of random DMIs")),
                        div(id="txtbox",textInput(inputId="text",label = "X-axis title", value = "Number of random DMIs")),
                        tags$style(type="text/css", "#txtbox {display: inline-block; max-width: 200px; }"),
                        div(id="txtbox", textInput("text2", label = "Y-axis title", value = "Frequency of random DMIs")),
                        div(id="txtbox",numericInput("xlimstart", label = "X-axis Start",0)),
                        div(id="txtbox",numericInput("xlimend", label = "X-axis End",200)),
                        tags$hr(),
                        tags$h4(tags$strong("Select Colors")),
                        
                        div(id="txtbox",colourpicker::colourInput("col", "Select bar colour", "firebrick3")),
                        div(id="txtbox",colourpicker::colourInput("col2", "Select background colour", "white")),
                        tags$hr(),
                        tags$h4(tags$strong("Select width/height to download plot as png")),
                        
                        div(id="txtbox",numericInput("width", label = "Width ", value = 1200)),
                        div(id="txtbox",numericInput("height", label = "Height ", value = 700))
                        
                    )),
                  tabPanel("Distribution of ELMs",
                           DT::dataTableOutput("diselmsdata"), tags$br(),tags$hr(),div(id="txtbox",actionButton("godis", "Interactive View"))
                  ),
                  tabPanel("Distribution of Domains",
                           DT::dataTableOutput("disdomdata"), tags$br(),tags$hr(),div(id="txtbox",actionButton("godisd", "Interactive View"))
                  ),
                  tabPanel("Network",fluidPage(tags$br(), selectInput("selectlayout", label = "Select Layout",
                                                                      choices = list("Circle" = "layout_in_circle","Nice" = "layout_nicely", "Random" = "layout_randomly", "Piecewise" = "piecewise.layout", "Gem" = "layout.gem"),
                                                                      selected = "piecewise.layout"),
                                               
                                               
                                               hr(),
                                               
                                               visNetworkOutput(outputId = "network",
                                                                height = "1500px",
                                                                width = "1500px")
                                               
                                               
                  )
                  )
                  
                  
      )
    )
  )),
  
  
  
  tabPanel(HTML("</a></li><li><a href=\"https://github.com/slimsuite/SLiMEnrich/wiki/Quick-Tutorial\", target = _blank>Getting Started")),tabPanel(HTML("</a></li><li><a href=\"https://github.com/slimsuite/SLiMEnrich/wiki\", target = _blank>Instructions"
  )), useShinyjs(),theme = shinytheme("sandstone"),
  tags$style(type="text/css", "#title {font-family: 'Impact', cursive;
             font-size: 32px;
             font-style:italic;
             font-color: #fff;
             -webkit-text-stroke-width: 1px;
             -webkit-text-stroke-color: black;}")
  
  
  ))

##############################
# Server logic
##############################
options(shiny.maxRequestSize=10000*1024^2)
server <- shinyServer(function(input, output, session){
  # This function computes a new data set. It can optionally take a function,
  # updateProgress, which will be called as each row of data is added.
  compute_data <- function(updateProgress = NULL) {
    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))
    
    for (i in 1:10) {
      Sys.sleep(0.25)
      
      # Compute new row of data
      new_row <- data.frame(x = rnorm(1), y = rnorm(1))
      
      # If we were passed a progress update function, call it
      if (is.function(updateProgress)) {
        text <- paste0("x:", round(new_row$x, 2), " y:", round(new_row$y, 2))
        updateProgress(detail = text)
      }
      
      # Add the new row of data
      dat <- rbind(dat, new_row)
    }
    
    dat
  }
  observeEvent(input$setting, {
    toggle(id = "settings", anim = TRUE)
  })
  observe({
    toggle(id = "settings")
  })
  observeEvent(input$uploadmotifs, {
    toggle(id = "uploadmotif", anim = TRUE)
  })
  observeEvent(input$SLiMrunid, {
    toggle(id = "slimrun", anim = TRUE)
    if(input$SLiMrunid){
      hide("slimf")
    }
    else{
      show("slimf")
    }
  })
  
  
  #lowercase function
  lowername <- function(x) {
    colnames(x) <- tolower(colnames(x))
    x
    
  }
  observeEvent(input$run, {
    MotifFile<-input$Motif
    PPIFile<-input$PPI
    if(is.null(PPIFile)){
      showNotification("PPI file is missing. Loading Example dataset", type = "error", duration = 5)
    }
    
    SliMJobId <- input$SLiMRun
    if(is.null(MotifFile) && SliMJobId == "" ){
      showNotification("SLiM file is missing. Loading Example Dataset", type = "error", duration = 5)
      showNotification("Loaded Example Dataset", type = "warning", duration = NULL)
      
    }
  })
  #####################################################Domain-Motif Interactions####################################################
  #*********************************************************************************************************************************
  #Uploaded Data
  ####################################################
  
  inputDataMotif <-eventReactive(input$run, {
    if(input$SLiMrunid){
      Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
    }
    else{
      MotifFile<-input$Motif
      
      if(is.null(MotifFile)){
        if(input$options == "pred"){
          #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
          Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
        }
        else if(input$options == "true") {
          
          Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
        }
      }
      
      else{
        #Read uploaded files
        Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      }
    }
  })
  inputDataPPI <-eventReactive(input$run, {
    PPIFile<-input$PPI
    if(is.null(PPIFile)){
      #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
      PPI2<-read.csv("example/adenofamilyPPIs.csv",header=TRUE,sep=",")
    }
    else{
      
      PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    }
    
  })
  
  inputDatadomain <-eventReactive(input$run, {
    #File upload check
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
      #dProtein <- lowername(dProtein)
    }
    
    else{
      dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
      # dProtein <- lowername(dProtein)
    }
  })
  inputDataMotifDomain <-eventReactive(input$run, {
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
      #Domain <- lowername(Domain)
      
    }
    else{
      Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
      #Domain <- lowername(Domain)
    }
    
  })
  
  #shows the data table
  output$udata<-renderDataTable({
    inputDataMotif()
    
    
  },
  caption = tags$h4(tags$strong("SLiM File"))
  )
  output$udata2<-renderDataTable({
    
    inputDataPPI()
    
  },
  caption = tags$h4(tags$strong("Interaction File"))
  )
  output$udata3<-renderDataTable({
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      return(NULL)
    }
    else{
      inputDatadomain()
    }
  },
  caption = tags$h4(tags$strong("Domain File"))
  )
  output$udata4<-renderDataTable({
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      return(NULL)
    }
    else{
      inputDataMotifDomain()
    }
    
  },
  caption = tags$h4(tags$strong("Motif-Domain File"))
  )
  #*************************************************************************************
  #Step 1: Potential DMIs
  ####################################################
  potentialDMIs <-eventReactive(input$run, {
    #File upload check
    MotifFile<-input$Motif
    
    
    #File upload check
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    
    else{
      dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
      
    }
    else{
      Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
    }
    if(input$SLiMrunid){
      Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
    }
    else{
      MotifFile<-input$Motif
      
      if(is.null(MotifFile)){
        if(input$options == "pred"){
          #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
          Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
        }
        else if(input$options == "true") {
          
          Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
        }
      }
      
      else{
        #Read uploaded files
        Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      }
    }
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    
    
    #Motif-Domain Mapping
    #Rename the columns in two files
    names(Motif_NR) <- c("Seq", "Motif")
    names(Domain) <- c("Motif", "Domain")
    if(input$withPro == TRUE){
      names(dProtein) <- c("Domain", "dProtein")
      #merge motif-domain and domain files
      elmPro <- merge(Domain, dProtein, by = "Domain")
      elmPro <- elmPro[, c("Motif","dProtein")]
      #Join/Merge two files based on Motif
      DMI <- merge(Motif_NR, elmPro, by="Motif")
      head(DMI)
      Uni_DMI <- unique(DMI)
      names(Uni_DMI) <- c("mProtein","Motif", "dProtein")
      Uni_DMI <- Uni_DMI[, c("mProtein" ,"Motif", "dProtein")]
      print(Uni_DMI)
    } else{
      
      #Join/Merge two files based on Motif
      Join <- merge(Motif_NR, Domain, by="Motif")
      names(Join) <- c("Motif", "Seq", "Domain")  
      #Domain-dProtein Mapping
      #Load results from the previous code)
      names(dProtein) <- c("Domain", "dProteins")
      #joined both files based on domain
      DMI <- merge(Join, dProtein,by="Domain")
      #Filtered unique DMIs
      Uni_DMI <- unique(DMI)
      #Named the header of output file
      names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
      Uni_DMI <- Uni_DMI[, c("mProtein","Motif", "Domain", "dProtein")]
      print(Uni_DMI)
    }
    return(Uni_DMI)
  })
  
  
  #shows the data table
  output$data<-DT::renderDataTable({
    #Run it only when run button is active
    if(input$run)
    {
      #Generates progress bar
      style <- isolate(input$style)
      
      # Create a Progress object
      progress <- shiny::Progress$new(style = style)
      progress$set(message = "Generating Data", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = "(Potential DMIs)")
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      formatStyle(datatable(potentialDMIs()), columns = 1:4, color = "black")
    }
    
  })
  #*************************************************************************************
  
  #Step 2: Predicted DMIs
  ####################################################
  predictedDMIs <-eventReactive(input$run, {
    if(input$SLiMrunid){
      Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
    }
    else{
      MotifFile<-input$Motif
      
      if(is.null(MotifFile)){
        if(input$options == "pred"){
          #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
          Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
        }
        else if(input$options == "true") {
          
          Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
        }
      }
      
      else{
        #Read uploaded files
        Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      }
    }
    PPIFile<-input$PPI
    if(is.null(PPIFile)){
      PPI2<-read.csv("example/adenofamilyPPIs.csv",header=TRUE,sep=",")
    }
    
    else{
      PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    }
    #File upload check
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    
    else{
      dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
      
    }
    else{
      Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
    }
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    
    #Rename the columns in two files
    names(Motif_NR) <- c("mProtein", "Motif")
    names(Domain) <- c("Motif", "Domain")
    
    if(input$withPro == TRUE){
      names(dProtein) <- c("Domain", "dProtein")
      #merge motif-domain and domain files
      elmPro <- merge(Domain, dProtein, by = "Domain")
      elmPro <- elmPro[, c("Motif","dProtein")]
      #Join/Merge two files based on Motif
      DMI <- merge(Motif_NR, elmPro, by="Motif")
      head(DMI)
      Uni_DMI <- unique(DMI)
      #names(Uni_DMI) <- c("mProtein","Motif", "dProtein")
      #Uni_DMI <- Uni_DMI[, c("mProtein" ,"Motif", "dProtein")]
      print(Uni_DMI)
      #PPI-DMI Mapping
      ########################################################################
      names(PPI2) <- c("mProtein", "dProtein")
      predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
      Uni_predDMIs <- unique(predDMI)
      #names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
      predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "dProtein")]
      print(predDMIs)
      
    } else{
      #Join/Merge two files based on Motif
      Join <- merge( Motif_NR, Domain, by="Motif")
      #print(Join)
      names(Join) <- c("Motif", "Seq", "Domains")  #Change header of the output file
      #Load mProtein_Motif_Domain file (result file from the previous code)
      names(dProtein) <- c("Domains", "dProteins")
      #joined both files based on Domains
      DMI <- merge(Join, dProtein,by="Domains")
      #Filtered unique DMIs
      Uni_DMI <- unique(DMI)
      #Named the header of output file
      names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
      #print(Uni_DMI)
      
      #PPI-DMI Mapping
      ########################################################################
      names(PPI2) <- c("mProtein", "dProtein")
      predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
      Uni_predDMIs <- unique(predDMI)
      names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
      predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
      print(predDMIs)
    }
    return(predDMIs)
  })
  #Shows predicted DMIs in DataTable
  output$PredDMIs<-DT::renderDataTable({
    #Run only if Run button is active
    if(input$run){
      #Progress bar
      style <- isolate(input$style)
      
      # Create a Progress object
      progress <- shiny::Progress$new(style = style)
      progress$set(message = "Generating Data", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = "(Predicted DMIs)")
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      
      
      formatStyle(datatable(predictedDMIs()), columns = 1:4, color = "black")
    }
  })
  #*************************************************************************************
  
  #Step 3: Statistics
  ####################################################
  summaryStat <- eventReactive(input$run, {
    if(input$SLiMrunid){
      Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
    }
    else{
      MotifFile<-input$Motif
      
      if(is.null(MotifFile)){
        if(input$options == "pred"){
          #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
          Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
        }
        else if(input$options == "true") {
          
          Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
        }
      }
      
      else{
        #Read uploaded files
        Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      }
    }
    PPIFile<-input$PPI
    if(is.null(PPIFile)){
      PPI2<-read.csv("example/adenofamilyPPIs.csv",header=TRUE,sep=",")
    }
    
    else{
      PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    }
    
    #File upload check
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    
    else{
      dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
      
    }
    else{
      Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
    }
    
    
    
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    #Rename the columns in two files
    #Rename the columns in two files
    names(Motif_NR) <- c("mProtein", "Motif")
    names(Domain) <- c("Motif", "Domain")
    
    if(input$withPro == TRUE){
      names(dProtein) <- c("Domain", "dProtein")
      #merge motif-domain and domain files
      elmPro <- merge(Domain, dProtein, by = "Domain")
      elmPro <- elmPro[, c("Motif","dProtein")]
      #Join/Merge two files based on Motif
      DMI <- merge(Motif_NR, elmPro, by="Motif")
      head(DMI)
      Uni_DMI <- unique(DMI)
      #names(Uni_DMI) <- c("mProtein","Motif", "dProtein")
      #Uni_DMI <- Uni_DMI[, c("mProtein" ,"Motif", "dProtein")]
      print(Uni_DMI)
      #PPI-DMI Mapping
      ########################################################################
      names(PPI2) <- c("mProtein", "dProtein")
      predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
      Uni_predDMIs <- unique(predDMI)
      #names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
      predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "dProtein")]
      print(predDMIs)
      colors=c("cadetblue1", "deepskyblue2", "darkblue")
      #Select unique Motif
      uniq_Motif <- unique(predDMI$Motif)
      a <- length(uniq_Motif)
      #Select unique mProtein
      uniq_mProtein <- unique(predDMI$mProtein)
      c <- length(uniq_mProtein)
      #Select unique dProteins
      uniq_dProtein <- unique(predDMI$dProtein)
      d <- length(uniq_dProtein)
      uniq_count<-c(a = length(uniq_Motif), c = length(uniq_mProtein), d = length(uniq_dProtein))
      #Create pie chart
      x <- c(a,c,d)
      
      statdmi <- data.frame(
        Datatype = factor(c("Motif","mProtein","dProtein")),
        Numbers = c(a,c,d)
      )
      
      p <- ggplot(data=statdmi, aes(x=Datatype, y=Numbers,fill=Datatype)) +
        geom_bar(colour="black", stat="identity") +
        guides(fill=FALSE)+
        scale_fill_manual(values = c("#85C1E9", "gold", "red"))
      
      p <- ggplotly(p)
      
    }else{
      #Join/Merge two files based on Motif
      Join <- merge( Motif_NR, Domain, by="Motif")
      #print(Join)
      names(Join) <- c("Motif", "Seq", "Domains")  #Change header of the output file
      #Load mProtein_Motif_Domain file (result file from the previous code)
      names(dProtein) <- c("Domains", "dProteins")
      #joined both files based on Domains
      DMI <- merge(Join, dProtein,by="Domains")
      #Filtered unique DMIs
      Uni_DMI <- unique(DMI)
      #Named the header of output file
      names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
      #print(Uni_DMI)
      #PPI-DMI Mapping
      ########################################################################
      names(PPI2) <- c("mProtein", "dProtein")
      predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
      Uni_predDMIs <- unique(predDMI)
      names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
      predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
      #print(predDMIs)
      colors=c("cadetblue1", "deepskyblue2", "blue", "darkblue")
      #Select unique Motif
      uniq_Motif <- unique(predDMI$Motif)
      a <- length(uniq_Motif)
      #Select unique Domain
      uniq_Domain <- unique(predDMI$Domain)
      b <- length(uniq_Domain)
      #Select unique mProtein
      uniq_mProtein <- unique(predDMI$mProtein)
      c <- length(uniq_mProtein)
      #Select unique dProteins
      uniq_dProtein <- unique(predDMI$dProtein)
      d <- length(uniq_dProtein)
      uniq_count<-c(a = length(uniq_Motif), b = length(uniq_Domain), c = length(uniq_mProtein), d = length(uniq_dProtein))
      #Create pie chart
      x <- c(a,b,c,d)
      
      statdmi <- data.frame(
        Datatype = factor(c("Motif","Domain","mProtein","dProtein")),
        Numbers = c(a,b,c,d)
      )
      
      p <- ggplot(data=statdmi, aes(x=Datatype, y=Numbers,fill=Datatype)) +
        geom_bar(colour="black", stat="identity") +
        guides(fill=FALSE)+
        scale_fill_manual(values = c("#9B59B6", "#85C1E9", "gold", "red"))
      
      p <- ggplotly(p)
      #p
      
      #barplot(x, main="Statistics of DMIs", col = colors)
      
      #legend("topright",
      #legend = c(paste("Motif=",a),paste("Domain=",b),paste("Motif containing Proteins=",c),paste("Domain containing Proteins=",d)), fill = colors)
    }
    p
  })
  output$plotbar <- renderPlotly({
    if(input$run){
      style <- isolate(input$style)
      
      # Create a Progress object
      progress <- shiny::Progress$new(style = style)
      progress$set(message ="Creating plot", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = "(Summary bar chart)")
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      summaryStat()
    }
  })
  #Step 4: Distribution of ELMs in the Predicted DMIs
  #####################################################
  #*************************************************************************************
  disELMs <- eventReactive(input$run, {
    if(input$SLiMrunid){
      Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
    }
    else{
      MotifFile<-input$Motif
      
      if(is.null(MotifFile)){
        if(input$options == "pred"){
          #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
          Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
        }
        else if(input$options == "true") {
          
          Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
        }
      }
      
      else{
        #Read uploaded files
        Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      }
    }
    PPIFile<-input$PPI
    if(is.null(PPIFile)){
      PPI2<-read.csv("example/adenofamilyPPIs.csv",header=TRUE,sep=",")
    }
    
    else{
      PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    }
    #File upload check
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    
    else{
      dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
      
    }
    else{
      Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
    }
    
    #Read uploaded files
    GOterms <- read.csv("data/elm_goterms.tsv",header=TRUE,sep="\t")
    names(GOterms) <- c("ELM", "GO Term", "Biological Function")
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    #Rename the columns in two files
    names(Motif_NR) <- c("mProtein", "Motif")
    names(Domain) <- c("Motif", "Domain")
    
    if(input$withPro == TRUE){
      names(dProtein) <- c("Domain", "dProtein")
      #merge motif-domain and domain files
      elmPro <- merge(Domain, dProtein, by = "Domain")
      elmPro <- elmPro[, c("Motif","dProtein")]
      #Join/Merge two files based on Motif
      DMI <- merge(Motif_NR, elmPro, by="Motif")
      head(DMI)
      Uni_DMI <- unique(DMI)
      #names(Uni_DMI) <- c("mProtein","Motif", "dProtein")
      #Uni_DMI <- Uni_DMI[, c("mProtein" ,"Motif", "dProtein")]
      print(Uni_DMI)
      #PPI-DMI Mapping
      ########################################################################
      names(PPI2) <- c("mProtein", "dProtein")
      predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
      Uni_predDMIs <- unique(predDMI)
      #names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
      predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "dProtein")]
      print(predDMIs)
      df_pred <- data.frame(predDMIs)["Motif"]
      names(df_pred) <- "Motif"
      print(df_pred)
      for (i in 1:length(df_pred)) {
        Matches <- count(df_pred)
        #names(Matches) <- "Frequency"
        #col <- cbind(df_pred,newcol)
        print(Matches)
        #
        
      }
      
      df_pred2 <- data.frame(Matches)[,(1:2)]
      names(df_pred2) <- c("ELM","Freq")
      print(df_pred2)
      Frequency <- df_pred2$Freq
      ELMs_names <- df_pred2$ELM
      df_pred2 <- data.frame(ELMs_names,Frequency)
      pvalueelm <- round(df_pred2$Frequency/nrow(df_pred),2)
      pvaluecol <- cbind(df_pred2,pvalueelm)
      names(pvaluecol) <- c("ELM", "Frequency", "Pvalue")
      GeneOntology <- merge(pvaluecol,GOterms, by="ELM")
      GeneOntology
      
    } else{
      #Join/Merge two files based on Motif
      Join <- merge( Motif_NR, Domain, by="Motif")
      #print(Join)
      names(Join) <- c("Motif", "Seq", "Domains")  #Change header of the output file
      #Load mProtein_Motif_Domain file (result file from the previous code)
      names(dProtein) <- c("Domains", "dProteins")
      #joined both files based on Domains
      DMI <- merge(Join, dProtein,by="Domains")
      #Filtered unique DMIs
      Uni_DMI <- unique(DMI)
      #Named the header of output file
      names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
      #PPI-DMI Mapping
      ########################################################################
      names(PPI2) <- c("mProtein", "dProtein")
      predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
      Uni_predDMIs <- unique(predDMI)
      names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
      #predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
      df_pred <- data.frame(Uni_predDMIs)["Motif"]
      names(df_pred) <- "Motif"
      print(df_pred)
      for (i in 1:length(df_pred)) {
        Matches <- count(df_pred)
        #names(Matches) <- "Frequency"
        #col <- cbind(df_pred,newcol)
        print(Matches)
        #
        
      }
      
      df_pred2 <- data.frame(Matches)[,(1:2)]
      names(df_pred2) <- c("ELM","Freq")
      print(df_pred2)
      Frequency <- df_pred2$Freq
      ELMs_names <- df_pred2$ELM
      df_pred2 <- data.frame(ELMs_names,Frequency)
      pvalueelm <- round(df_pred2$Frequency/nrow(df_pred),2)
      pvaluecol <- cbind(df_pred2,pvalueelm)
      names(pvaluecol) <- c("ELM", "Frequency", "Pvalue")
      GeneOntology <- merge(pvaluecol,GOterms, by="ELM")
      GeneOntology
    }
  })
  
  output$diselmsdata <-DT::renderDataTable({
    #Run only if Run button is active
    if(input$run){
      #Progress bar
      style <- isolate(input$style)
      
      # Create a Progress object
      progress <- shiny::Progress$new(style = style)
      progress$set(message = "Generating Data", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = "ELM distribution")
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      
      disELMs()
      
    }
  })
  displotfunc <- function(){
    if(input$SLiMrunid){
      Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
    }
    else{
      MotifFile<-input$Motif
      
      if(is.null(MotifFile)){
        if(input$options == "pred"){
          #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
          Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
        }
        else if(input$options == "true") {
          
          Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
        }
      }
      
      else{
        #Read uploaded files
        Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      }
    }
    PPIFile<-input$PPI
    if(is.null(PPIFile)){
      PPI2<-read.csv("example/adenofamilyPPIs.csv",header=TRUE,sep=",")
    }
    
    else{
      PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    }
    #File upload check
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    
    else{
      dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
      
    }
    else{
      Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
    }
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    #Rename the columns in two files
    names( Motif_NR) <- c("Seq", "Motif")
    names(Domain) <- c("Motif", "Domain")
    #Join/Merge two files based on Motif
    Join <- merge( Motif_NR, Domain, by="Motif")
    #print(Join)
    names(Join) <- c("Motif", "Seq", "Domains")  #Change header of the output file
    #Load mProtein_Motif_Domain file (result file from the previous code)
    names(dProtein) <- c("Domains", "dProteins")
    #joined both files based on Domains
    DMI <- merge(Join, dProtein,by="Domains")
    #Filtered unique DMIs
    Uni_DMI <- unique(DMI)
    #Named the header of output file
    names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
    #PPI-DMI Mapping
    ########################################################################
    names(PPI2) <- c("mProtein", "dProtein")
    predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
    Uni_predDMIs <- unique(predDMI)
    names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
    #predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
    df_pred <- data.frame(Uni_predDMIs)["Motif"]
    names(df_pred) <- "Motif"
    print(df_pred)
    for (i in 1:length(df_pred)) {
      Matches <- count(df_pred)
      #names(Matches) <- "Frequency"
      #col <- cbind(df_pred,newcol)
      print(Matches)
      #
      
    }
    df_pred2 <- data.frame(Matches)[,(1:2)]
    names(df_pred2) <- c("ELM","Freq")
    print(df_pred2)
    Frequency <- df_pred2$Freq
    ELMs_names <- df_pred2$ELM
    df_pred2 <- data.frame(ELMs_names,Frequency)
    disdmi <- data.frame(
      Datatype = factor(df_pred2$ELM),
      Numbers = df_pred2$Freq
    )
    #Frequency <- df_pred2$Freq
    #ELMs_names <- df_pred2$ELM
    #write.csv(df_pred2, "ELMs.csv", row.names = FALSE)
    # A basic box with the conditions colored
    #ggplot(df_pred2, aes(depth, fill = cut)) +
    # geom_density(position = "stack")
    
    p <- ggplot(data=disdmi, aes(x=Datatype, y=Numbers,fill=Datatype)) +
      geom_bar(colour="black", stat="identity") +
      guides(fill=FALSE)+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
    p <- ggplotly(p)
    
  }
  
  output$diselmchart <- renderPlotly({
    #Run only if Run button is active
    if(input$godis){
      #Progress bar
      style <- isolate(input$style)
      
      # Create a Progress object
      progress <- shiny::Progress$new(style = style)
      progress$set(message = "Generating Chart", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = "Generating interactive view")
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      displotfunc()
    }
  })
  #Distribution of Domains in the Predicted DMIs
  #####################################################
  #*************************************************************************************
  disDom <- eventReactive(input$run, {
    if(input$SLiMrunid){
      Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
    }
    else{
      MotifFile<-input$Motif
      
      if(is.null(MotifFile)){
        if(input$options == "pred"){
          #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
          Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
        }
        else if(input$options == "true") {
          
          Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
        }
      }
      
      else{
        #Read uploaded files
        Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      }
    }
    PPIFile<-input$PPI
    if(is.null(PPIFile)){
      PPI2<-read.csv("example/adenofamilyPPIs.csv",header=TRUE,sep=",")
    }
    
    else{
      PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    }
    #File upload check
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    
    else{
      dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
      
    }
    else{
      Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
    }
    
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    #Rename the columns in two files
    names( Motif_NR) <- c("Seq", "Motif")
    names(Domain) <- c("Motif", "Domain")
    #Join/Merge two files based on Motif
    Join <- merge( Motif_NR, Domain, by="Motif")
    #print(Join)
    names(Join) <- c("Motif", "Seq", "Domains")  #Change header of the output file
    #Load mProtein_Motif_Domain file (result file from the previous code)
    names(dProtein) <- c("Domains", "dProteins")
    #joined both files based on Domains
    DMI <- merge(Join, dProtein,by="Domains")
    #Filtered unique DMIs
    Uni_DMI <- unique(DMI)
    #Named the header of output file
    names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
    #PPI-DMI Mapping
    ########################################################################
    names(PPI2) <- c("mProtein", "dProtein")
    predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
    Uni_predDMIs <- unique(predDMI)
    names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
    #predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
    df_pred <- data.frame(Uni_predDMIs)["Domain"]
    names(df_pred) <- "Domain"
    print(df_pred)
    for (i in 1:length(df_pred)) {
      Matches <- count(df_pred)
      #names(Matches) <- "Frequency"
      #col <- cbind(df_pred,newcol)
      print(Matches)
      #
      
    }
    
    df_pred2 <- data.frame(Matches)[,(1:2)]
    names(df_pred2) <- c("Domain","Freq")
    print(df_pred2)
    Frequency <- df_pred2$Freq
    ELMs_names <- df_pred2$Domain
    df_pred2 <- data.frame(ELMs_names,Frequency)
    pvalueelm <- round(df_pred2$Frequency/nrow(df_pred),2)
    pvaluecol <- cbind(df_pred2,pvalueelm)
    names(pvaluecol) <- c("Domain", "Frequency", "Pvalue")
    pvaluecol
  })
  
  output$disdomdata <-DT::renderDataTable({
    #Run only if Run button is active
    if(input$run){
      #Progress bar
      style <- isolate(input$style)
      
      # Create a Progress object
      progress <- shiny::Progress$new(style = style)
      progress$set(message = "Generating Data", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = "Domain distribution")
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      
      disDom()
      
    }
  })
  disdomplotfunc <- function(){
    if(input$SLiMrunid){
      Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
    }
    else{
      MotifFile<-input$Motif
      
      if(is.null(MotifFile)){
        if(input$options == "pred"){
          #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
          Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
        }
        else if(input$options == "true") {
          
          Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
        }
      }
      
      else{
        #Read uploaded files
        Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      }
    }
    PPIFile<-input$PPI
    if(is.null(PPIFile)){
      PPI2<-read.csv("example/adenofamilyPPIs.csv",header=TRUE,sep=",")
    }
    
    else{
      PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    }
    #File upload check
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    
    else{
      dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
      
    }
    else{
      Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
    }
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    #Rename the columns in two files
    names( Motif_NR) <- c("Seq", "Motif")
    names(Domain) <- c("Motif", "Domain")
    #Join/Merge two files based on Motif
    Join <- merge( Motif_NR, Domain, by="Motif")
    #print(Join)
    names(Join) <- c("Motif", "Seq", "Domains")  #Change header of the output file
    #Load mProtein_Motif_Domain file (result file from the previous code)
    names(dProtein) <- c("Domains", "dProteins")
    #joined both files based on Domains
    DMI <- merge(Join, dProtein,by="Domains")
    #Filtered unique DMIs
    Uni_DMI <- unique(DMI)
    #Named the header of output file
    names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
    #PPI-DMI Mapping
    ########################################################################
    names(PPI2) <- c("mProtein", "dProtein")
    predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
    Uni_predDMIs <- unique(predDMI)
    names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
    #predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
    df_pred <- data.frame(Uni_predDMIs)["Domain"]
    names(df_pred) <- "Domain"
    print(df_pred)
    for (i in 1:length(df_pred)) {
      Matches <- count(df_pred)
      #names(Matches) <- "Frequency"
      #col <- cbind(df_pred,newcol)
      print(Matches)
      #
      
    }
    df_pred2 <- data.frame(Matches)[,(1:2)]
    names(df_pred2) <- c("Domain","Freq")
    print(df_pred2)
    Frequency <- df_pred2$Freq
    dom_names <- df_pred2$Domain
    df_pred2 <- data.frame(dom_names,Frequency)
    disdomdmi <- data.frame(
      Datatype = factor(df_pred2$dom_names),
      Numbers = df_pred2$Frequency
    )
    
    p <- ggplot(data=disdomdmi, aes(x=Datatype, y=Numbers,fill=Datatype)) +
      geom_bar(colour="black", stat="identity") +
      guides(fill=FALSE)+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
    p <- ggplotly(p)
    
    
  }
  
  output$disdomchart <- renderPlotly({
    #Run only if Run button is active
    if(input$godisd){
      #Progress bar
      style <- isolate(input$style)
      
      # Create a Progress object
      progress <- shiny::Progress$new(style = style)
      progress$set(message = "Generating Chart", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = "Generating interactive view")
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      disdomplotfunc()
    }
  })
  
  #####################################################Enrichment Analysis####################################################
  #*********************************************************************************************************************************
  
  
  #*************************************************************************************
  
  #Step 5: rPPI-DMI Mapping
  ####################################################
  #This step generates 1000 random files and then compares each file with the potential DMIs dataset to generate a list of numbers (matches found in each PPI file).
  #Randomized files will be stored in a new directory created in App folder named as "RandomFiles" and can be accessed later if required.
  #A file named randomNumbers will be generated in "RandomFiles" that will be used to generate Histogram later.
  ####################################################
  rPPIDMI <-reactive({
    if(input$SLiMrunid){
      Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
    }
    else{
      MotifFile<-input$Motif
      
      if(is.null(MotifFile)){
        if(input$options == "pred"){
          #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
          Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
        }
        else if(input$options == "true") {
          
          Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
        }
      }
      
      else{
        #Read uploaded files
        Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      }
    }
    PPIFile<-input$PPI
    if(is.null(PPIFile)){
      PPI<-read.csv("example/adenofamilyPPIs.csv",header=TRUE,sep=",")
    }
    
    else{
      PPI<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    }
    #File upload check
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    
    else{
      dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
      
    }
    else{
      Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
    }
    observe({
      show(id = "go2")
    })
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    #Rename the columns in two files
    names(Motif_NR) <- c("mProtein", "Motif")
    names(Domain) <- c("Motif", "Domain")
    
    if(input$withPro == TRUE){
      names(dProtein) <- c("Domain", "dProtein")
      #merge motif-domain and domain files
      elmPro <- merge(Domain, dProtein, by = "Domain")
      elmPro <- elmPro[, c("Motif","dProtein")]
      #Join/Merge two files based on Motif
      DMI <- merge(Motif_NR, elmPro, by="Motif")
      head(DMI)
      Uni_DMI <- unique(DMI)
      #names(Uni_DMI) <- c("mProtein","Motif", "dProtein")
      #Uni_DMI <- Uni_DMI[, c("mProtein" ,"Motif", "dProtein")]
      print(Uni_DMI)
      #PPI-DMI Mapping
      ########################################################################
      names(PPI) <- c("mProtein", "dProtein")
      predDMI <- merge(PPI, Uni_DMI, by= c("mProtein", "dProtein"))
      Uni_predDMIs <- unique(predDMI)
      #names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
      predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "dProtein")]
      print(predDMIs)
      #Randomization/Permutations
      ##############################################################################
      PPI_data <- PPI
      PPI_data <- unique(PPI)
      names(PPI_data) <- c("mProtein", "dProtein")
      PPI_Matrix<-matrix(data = PPI_data$mProtein)
      PPI_Matrix2<-matrix(data = PPI_data$dProtein)
      PermutationFunction <- function (data, k) {
        
        # creating matrix: amount of variables * amount of permutations
        permutations <- matrix(1:(k * length(data[1, ])), nrow=k)
        row <- NULL
        
        # For loop runs for i number of times from 1:length(data[1,])) (*which is number of rows in the dataset*).
        #sample(data[ , i] randomizes the column
        #k is the number of entries in the dataset
        for (i in 1:length(data[1,])) {
          permutations[ ,i] <- sample(data[ , i], k, replace=FALSE)
        }
        permutations
      }
      #dir.create("RandomFiles")
      
      #dirName <- paste0("RandomFiles_", strsplit(as.character(PPIFile$name), '.csv'))
      #dir.create(dirName)
      #for loop to create 1000 randomized files
      showNotification("Performing the randomizations", type = "message", duration = 5)
      data <- list()
      for (j in 1:1000) {
        permutation<-PermutationFunction(PPI_Matrix, k = length(PPI_Matrix))
        permutation2<-PermutationFunction(PPI_Matrix2, k = length(PPI_Matrix2))
        final_file<- c(paste(permutation,permutation2, sep = ":"))
        newCol1<-strsplit(as.character(final_file),':',fixed=TRUE)
        df<-data.frame(final_file,do.call(rbind, newCol1))
        subset = df[,c(2,3)]
        names(subset)<-c("mProtein", "dProtein")
        data[[j]] <- subset
        #head(data[[j]])
        #write.csv(perm, "rPPI.csv", row.names = FALSE)
      }
      #rPPI-DMI Mapping                                                               
      #################################################################################
      showNotification("1000 random PPI files have been created", type = "message", duration = 5)
      showNotification("Now predicing DMIs from the random PPI data", type = "message", closeButton = TRUE,duration = 15)
      m <- data.frame()
      for (i in 1:1000) {
        rPPI <- data[[i]]
        names(rPPI)<-c("mProtein", "dProtein")
        #names(Uni_DMI) <- c("Domain", "Motif", "MotifProtein", "DomainProtein")
        DMI_rPPI <- merge(Uni_DMI, rPPI, by= c("mProtein", "dProtein"))
        Matches <- nrow(DMI_rPPI)
        print(Matches)
        dmatch <- data.frame(Matches)
        m=rbind(m,dmatch)
        row.names(m) <- NULL
      }
      
    } else{
      
      #Join/Merge two files based on Motif
      Join <- merge( Motif_NR, Domain, by="Motif")
      #print(Join)
      names(Join) <- c("Motif", "Seq", "Domains")  #Change header of the output file
      #Load mProtein_Motif_Domain file (result file from the previous code)
      names(dProtein) <- c("Domains", "dProteins")
      #joined both files based on Domains
      DMI <- merge(Join, dProtein,by="Domains")
      #Filtered unique DMIs
      Uni_DMI <- unique(DMI)
      #Named the header of output file
      names(PPI) <- c("mProtein", "dProtein")
      names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
      predDMI <- merge(PPI, Uni_DMI, by= c("mProtein", "dProtein"))
      Uni_predDMIs <- unique(predDMI)
      names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
      predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
      #print(predDMIs)
      #Randomization/Permutations
      ##############################################################################
      PPI_data <- PPI
      PPI_data <- unique(PPI)
      names(PPI_data) <- c("mProtein", "dProtein")
      PPI_Matrix<-matrix(data = PPI_data$mProtein)
      PPI_Matrix2<-matrix(data = PPI_data$dProtein)
      PermutationFunction <- function (data, k) {
        
        # creating matrix: amount of variables * amount of permutations
        permutations <- matrix(1:(k * length(data[1, ])), nrow=k)
        row <- NULL
        
        # For loop runs for i number of times from 1:length(data[1,])) (*which is number of rows in the dataset*).
        #sample(data[ , i] randomizes the column
        #k is the number of entries in the dataset
        for (i in 1:length(data[1,])) {
          permutations[ ,i] <- sample(data[ , i], k, replace=FALSE)
        }
        permutations
      }
      #dir.create("RandomFiles")
      
      #dirName <- paste0("RandomFiles_", strsplit(as.character(PPIFile$name), '.csv'))
      #dir.create(dirName)
      #for loop to create 1000 randomized files
      showNotification("Performing the randomizations", type = "message", duration = 5)
      data <- list()
      for (j in 1:1000) {
        permutation<-PermutationFunction(PPI_Matrix, k = length(PPI_Matrix))
        permutation2<-PermutationFunction(PPI_Matrix2, k = length(PPI_Matrix2))
        final_file<- c(paste(permutation,permutation2, sep = ":"))
        newCol1<-strsplit(as.character(final_file),':',fixed=TRUE)
        df<-data.frame(final_file,do.call(rbind, newCol1))
        subset = df[,c(2,3)]
        names(subset)<-c("mProtein", "dProtein")
        data[[j]] <- subset
        #head(data[[j]])
        #write.csv(perm, "rPPI.csv", row.names = FALSE)
      }
      #rPPI-DMI Mapping                                                               
      #################################################################################
      showNotification("1000 random PPI files have been created", type = "message", duration = 5)
      showNotification("Now predicing DMIs from the random PPI data", type = "message", closeButton = TRUE,duration = 15)
      m <- data.frame()
      for (i in 1:1000) {
        rPPI <- data[[i]]
        names(rPPI)<-c("MotifProtein", "DomainProtein")
        names(Uni_DMI) <- c("Domain", "Motif", "MotifProtein", "DomainProtein")
        DMI_rPPI <- merge(Uni_DMI, rPPI, by= c("MotifProtein", "DomainProtein"))
        Matches <- nrow(DMI_rPPI)
        print(Matches)
        dmatch <- data.frame(Matches)
        m=rbind(m,dmatch)
        row.names(m) <- NULL
      }
    }
    #output_randomNumbers <- write.csv(m, paste0("randomNumbers.csv"), row.names = FALSE)
    m
  })
  #*************************************************************************************
  
  #Step 6: Histogram
  ####################################################
  
  output$histogram <- renderPlot({
    
    if(input$run){
      style <- isolate(input$style)
      
      # Create a Progress object
      progress <- shiny::Progress$new(style = style)
      progress$set(message ="Creating Plot", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = "This may take a while!!. Nothing will respond while it's being calculated (e.g. network)")
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      rPPIDMI()
      plotInput()
    }
    
  })
  #Function to generate Histogram
  plotInput <- function(){
    if(input$SLiMrunid){
      Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
    }
    else{
      MotifFile<-input$Motif
      
      if(is.null(MotifFile)){
        if(input$options == "pred"){
          #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
          Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
        }
        else if(input$options == "true") {
          
          Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
        }
      }
      
      else{
        #Read uploaded files
        Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      }
    }
    PPIFile<-input$PPI
    if(is.null(PPIFile)){
      PPI<-read.csv("example/adenofamilyPPIs.csv",header=TRUE,sep=",")
    }
    
    else{
      PPI<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    }
    #File upload check
    DomainFile<-input$domain
    if(is.null(DomainFile)){
      dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    
    else{
      dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
      dProtein <- lowername(dProtein)
      dProtein <- dProtein[,c("pfam","accnum")]
    }
    MotifDomainFile<-input$MotifDomain
    if(is.null(MotifDomainFile)){
      Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
      
    }
    else{
      Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
      Domain <- lowername(Domain)
      Domain <- Domain[, c("elmidentifier","interactiondomainid")]
    }
    
    #dirName <- paste0("RandomFiles_", strsplit(as.character(PPIFile$name), '.csv'))
    x <- rPPIDMI()
    names(x) <- "values"
    bins<- seq(min(x), max(x), length.out = input$bins + 1)
    par(bg=input$col2)
    h<- hist(x$values, breaks=bins, col = input$col, border = 'black', main=input$text3, ylab=input$text2,
             xlab=input$text, xlim = c( input$xlimstart,  input$xlimend), labels = input$barlabel, cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
    xfit <- seq(min(x$values),max(x$values),length=40)
    
    yfit <- dnorm(xfit, mean=mean(x$values),sd=sd(x$values))
    
    
    yfit<- yfit*diff(h$mids[1:2])*length(x$values)
    
    
    lines(xfit,yfit)
    ob_fdr <- x[x$values <= nrow(predictedDMIs()), ]
    #axis(side=3, lwd = 0, lwd.ticks = 4, at=nrow(predictedDMIs()), lend=1, labels = FALSE, tcl=5, font=2, col = "black", padj = 0, lty = 3)
    #shows the observed value
    mtext(paste("Observed value is: ", nrow(predictedDMIs())), side = 3, at=nrow(predictedDMIs()), font = 4)
    mtext(paste0("P-value is: ", length(x[x >= nrow(predictedDMIs())])/1000), side = 3, at=nrow(predictedDMIs())+90, font = 4, col = "red")
    #mtext(paste0("Mean is: ", mean(x$values)), side = 3, at=mean(x$values), font = 4, col="red")
    #mtext(paste0("FDR is: ", mean(x$values)/nrow(predictedDMIs())), side = 3, at=50, font = 4, col= "red")
    #points arrow on the observed value
    arrows(nrow(predictedDMIs()), 480, nrow(predictedDMIs()), 0, lwd = 2, col = "black", length = 0.1, lty = 3)
    pvalue <-  paste0("<b>P-value is: </b>", length(x[x >= nrow(predictedDMIs())])/1000)
    meanvalue <- paste0("<b>Mean is: </b>", round(mean(x$values)))
    Escore <- paste0("<b>Enrichment score (E-score) is: </b>", round(nrow(predictedDMIs())/mean(x$values),2))
    FDR <- paste0("<b>False Discrovery Rate (FDR) is: </b>", round(mean(ob_fdr)/nrow(predictedDMIs()),2))
    output$summary <- renderUI({
      
      HTML(paste("<font color=\"#FF0000\"><b>Summary of Histogram</b></font>", pvalue, meanvalue, Escore, FDR, sep = '<hr/>'))
      
    })
  }
  compute_data <- function(updateProgress = NULL) {
    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))
    
    for (i in 1:10) {
      Sys.sleep(0.25)
      
      # Compute new row of data
      new_row <- data.frame(x = rnorm(1), y = rnorm(1))
      
      # If we were passed a progress update function, call it
      if (is.function(updateProgress)) {
        text <- paste0("x:", round(new_row$x, 2), " y:", round(new_row$y, 2))
        updateProgress(detail = text)
      }
      
      # Add the new row of data
      dat <- rbind(dat, new_row)
    }
    
    dat
  }
  
  #creates plot in  a seperate window
  output$plot <- renderPlot({
    plotInput()
  })
  
  #*************************************************************************************
  
  #Step 7: DMI Network
  ####################################################
  mynetwork <- function(){
    set.seed(20)
    if(input$run){
      style <- isolate(input$style)
      
      # Create a Progress object
      progress <- shiny::Progress$new(style = style)
      progress$set(message ="Creating Network", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a closure to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = "")
      }
      
      # Compute the new data, and pass in the updateProgress function so
      # that it can update the progress indicator.
      compute_data(updateProgress)
      if(input$SLiMrunid){
        Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
      }
      else{
        MotifFile<-input$Motif
        
        if(is.null(MotifFile)){
          if(input$options == "pred"){
            #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
            Motif<-read.csv("example/slimprob.occ.csv",header=TRUE,sep=",")
          }
          else if(input$options == "true") {
            
            Motif<-read.csv("data/known.occ.csv",header=TRUE,sep=",")
          }
        }
        
        else{
          #Read uploaded files
          Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
        }
      }
      PPIFile<-input$PPI
      if(is.null(PPIFile)){
        PPI2<-read.csv("example/adenofamilyPPIs.csv",header=TRUE,sep=",")
      }
      
      else{
        PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
      }
      #File upload check
      DomainFile<-input$domain
      if(is.null(DomainFile)){
        dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")
        dProtein <- lowername(dProtein)
        dProtein <- dProtein[,c("pfam","accnum")]
      }
      
      else{
        dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")
        dProtein <- lowername(dProtein)
        dProtein <- dProtein[,c("pfam","accnum")]
      }
      MotifDomainFile<-input$MotifDomain
      if(is.null(MotifDomainFile)){
        Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")
        Domain <- lowername(Domain)
        Domain <- Domain[, c("elmidentifier","interactiondomainid")]
        
      }
      else{
        Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")
        Domain <- lowername(Domain)
        Domain <- Domain[, c("elmidentifier","interactiondomainid")]
      }
      #Read uploaded files
      Motif <- lowername(Motif)
      Motif <- Motif[, c("accnum","motif")]
      names(Motif) <- c("UniprotID","Motif")
      Motif_NR<-unique(Motif)
      #Rename the columns in two files
      names(Motif_NR) <- c("mProtein", "Motif")
      names(Domain) <- c("Motif", "Domain")
      
      if(input$withPro == TRUE){
        names(dProtein) <- c("Domain", "dProtein")
        #merge motif-domain and domain files
        elmPro <- merge(Domain, dProtein, by = "Domain")
        elmPro <- elmPro[, c("Motif","dProtein")]
        #Join/Merge two files based on Motif
        DMI <- merge(Motif_NR, elmPro, by="Motif")
        head(DMI)
        Uni_DMI <- unique(DMI)
        #names(Uni_DMI) <- c("mProtein","Motif", "dProtein")
        #Uni_DMI <- Uni_DMI[, c("mProtein" ,"Motif", "dProtein")]
        print(Uni_DMI)
        #PPI-DMI Mapping
        ########################################################################
        names(PPI2) <- c("mProtein", "dProtein")
        predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
        Uni_predDMIs <- unique(predDMI)
        #names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
        predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "dProtein")]
        print(predDMIs)
        if(nrow(predDMIs) != 0){
          #Network
          
          first <- predDMI[,c("mProtein","Motif")]
          g <- graph.data.frame(first, directed = F)
          #igraph.options(plot.layout=layout.graphopt, vertex.size=10)
          V(g)$color <- ifelse(V(g)$name %in% predDMI[,1], "#A93226", "#F7DC6F")
          V(g)$shape <- ifelse(V(g)$name %in% predDMI[,1], "square", "box")
          
          second <- predDMI[,c("Motif","dProtein")]
          g2 <- graph.data.frame(second, directed = F)
          #igraph.options(plot.layout=layout.graphopt, vertex.size=10)
          V(g2)$color <- ifelse(V(g2)$name %in% predDMI[,2], "#85C1E9", "#85C1E9")
          V(g2)$shape <- ifelse(V(g2)$name %in% predDMI[,2], "circle", "vrectangle")
          
          
          #merge networks
          g4 = graph.union(g,g2, byname = TRUE)
          #g4 <- g %u% g2 %u% g3
          g5 <- simplify(g4, remove.multiple = TRUE, remove.loops = TRUE)
          #plot.igraph(g5, layout=layout.fruchterman.reingold)
          V(g5)$color <- ifelse(is.na(V(g5)$color_1),
                                V(g5)$color_2,V(g5)$color_1)
          V(g5)$shape <- ifelse(is.na(V(g5)$shape_1),
                                V(g5)$shape_2,V(g5)$shape_1)
          E(g5)$color <- "black"
          E(g)$width <- 9
          g6 <- visIgraph(g5,layout = input$selectlayout, physics = FALSE, smooth = TRUE, type = "square")
          #visExport(g6, type = "png", name = "export-network",float = "left", label = "Save network", background = "white", style= "")
          
        }
      }else{
        
        #Join/Merge two files based on Motif
        Join <- merge( Motif_NR, Domain, by="Motif")
        #print(Join)
        names(Join) <- c("Motif", "Seq", "Domains")  #Change header of the output file
        #Load mProtein_Motif_Domain file (result file from the previous code)
        names(dProtein) <- c("Domains", "dProteins")
        #joined both files based on Domains
        DMI <- merge(Join, dProtein,by="Domains")
        #Filtered unique DMIs
        Uni_DMI <- unique(DMI)
        #Named the header of output file
        names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
        #print(Uni_DMI)
        
        #PPI-DMI Mapping
        ########################################################################
        names(PPI2) <- c("mProtein", "dProtein")
        predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
        
        Uni_predDMIs <- unique(predDMI)
        names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
        predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
        if(nrow(predDMIs) != 0){
          #Network
          
          first <- predDMI[,c("mProtein","Motif")]
          g <- graph.data.frame(first, directed = F)
          #igraph.options(plot.layout=layout.graphopt, vertex.size=10)
          V(g)$color <- ifelse(V(g)$name %in% predDMI[,1], "#A93226", "#F7DC6F")
          V(g)$shape <- ifelse(V(g)$name %in% predDMI[,1], "square", "box")
          #plot(g,  edge.color="orange")
          #visIgraph(g)
          ##print(first)
          #V(g)$color <- "red"
          #Second
          second <- predDMI[,c("Motif","Domain")]
          g2 <- graph.data.frame(second, directed = F)
          #igraph.options(plot.layout=layout.graphopt, vertex.size=10)
          V(g2)$color <- ifelse(V(g2)$name %in% predDMI[,1], "#9B59B6", "#D35400")
          V(g2)$shape <- ifelse(V(g2)$name %in% predDMI[,1], "box", "circle")
          #plot(g2,  edge.color="orange")
          #print(second)
          #visIgraph(g2)
          #V(g2)$color <- "green"
          #Third
          
          third <- predDMI[,c("Domain","dProtein")]
          g3 <- graph.data.frame(third, directed = F)
          #igraph.options(plot.layout=layout.graphopt, vertex.size=10)
          V(g3)$color <- ifelse(V(g3)$name %in% predDMI[,3], "#9B59B6", "#85C1E9")
          V(g3)$shape <- ifelse(V(g3)$name %in% predDMI[,3], "vrectangle", "circle")
          #plot(g3,  edge.color="orange")
          #visIgraph(g3)
          #print(third)
          #V(g3)$color <- "pink"
          
          #merge networks
          g4 = graph.union(g,g2,g3, byname = TRUE)
          #g4 <- g %u% g2 %u% g3
          g5 <- simplify(g4, remove.multiple = TRUE, remove.loops = TRUE)
          #plot.igraph(g5, layout=layout.fruchterman.reingold)
          V(g5)$color <- ifelse(is.na(V(g5)$color_1),
                                V(g5)$color_3,V(g5)$color_1)
          V(g5)$shape <- ifelse(is.na(V(g5)$shape_1),
                                V(g5)$shape_3,V(g5)$shape_1)
          E(g5)$color <- "black"
          E(g)$width <- 9
          g6 <- visIgraph(g5,layout = input$selectlayout, physics = FALSE, smooth = TRUE, type = "square")
          #visExport(g6, type = "png", name = "export-network",float = "left", label = "Save network", background = "white", style= "")
          
        }
      }
    }
  }
  output$network <- renderVisNetwork({
    mynetwork()
    
  })
  
  ###
  #Buttons
  #**************
  #download button to download potential DMIs data
  output$downloadDMI <- downloadHandler(
    filename = "potentialDMIs.csv",
    content = function(file) {
      write.csv(potentialDMIs(), file)
    }
  )
  #download button to download predicted DMIs data
  output$downloadpredDMI <- downloadHandler(
    filename = "predDMIs.csv",
    content = function(file) {
      write.csv(predictedDMIs(), file)
    }
  )
  #download button to download the plot
  output$downloadPlot <- downloadHandler(
    filename = 'Histogram.png',
    content = function(file) {
      png(file, width = input$width, height = input$height, units = "px", pointsize = 12, res = 300)
      plotInput()
      dev.off()
    }
  )
  #download button to download the plot
  output$downloadrandom <- downloadHandler(
    filename = "predDMIs.csv",
    content = function(file) {
      write.csv(predictedDMIs(), file)
    }
  )
})


shinyApp(ui = ui, server = server)



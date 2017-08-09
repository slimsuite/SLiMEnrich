#*********************************************************************************************************
#*********************************************************************************************************
# Short Linear Motif Enrichment Analysis App (SLiMEnrich)
# Developer: **Sobia Idrees**
# Version: 1.0.7
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
#V1.0.6 - Added p-values, GOTerms and Biological function of ELMs in ELM distribution tab
#V1.0.7 - File headers to lowercase for consistency
##############################
#SLiMEnrich program is free software: you can redistribute it and/or modify
 #   it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation, either version 3 of the License, or
  #  (at your option) any later version.

   # SLiMEnrich program is distributed in the hope that it will be useful,
   # but WITHOUT ANY WARRANTY; without even the implied warranty of
   # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   # GNU General Public License for more details.

   # You should have received a copy of the GNU General Public License
   # along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################
#Required Libraries
##############################
# Check whether packages of interest are installed
is_installed = function(mypkg) is.element(mypkg, installed.packages()[,1]) 
# Install library if not already installed
# Run a for-loop of all the package names listed below in the function call
# with the list of packages: load_or_install(c("pkg1", "pkg2",..., "pkgn"))
load_or_install = function(package_names) 
{ 
  for(package_name in package_names) 
  { 
    if(!is_installed(package_name)) 
    { 
      #install.packages(package_name,repos="http://lib.stat.cmu.edu/R/CRAN") 
      install.packages(package_name)
    } 
    library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE) 
  } 
}
load_or_install(c("shiny", "ggplot2", "colourpicker", "shinyBS", "shinythemes", "DT", "shinyjs", "visNetwork", "igraph","markdown","plotly", "plyr"))
#library(shiny)
#library(ggplot2)
#library(colourpicker)
#library(shinyBS)
#library(shinythemes)
#library(DT)
#library(shinyjs)
#library(visNetwork)
#library(igraph)
#library(plotly)
#library(markdown)
##############################
#GUI of the App
##############################
if (interactive()) {
#navbar page with sidebar layout along with tabsets
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
                    "))
    ),
    # Sidebar
    sidebarLayout(
      sidebarPanel(
        fileInput("PPI","Select Interaction file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
        
        fileInput("Motif","Select SLiM prediction file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
        actionButton("run", "Run", width = "100px"),
        div(id="fileuploads",checkboxInput("uploadmotifs",label = "Upload Domain and/or Motif-Domain Files", value = FALSE)),
        div(id="uploadmotif",  fileInput("domain","Select Domain file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
            fileInput("MotifDomain","Select Motif-Domain file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')))
      ),
      
      # MainPanel
      mainPanel(
        #Creates a seperate window (pop up window)
        #Creates a seperate window (pop up window)
        bsModal("DisE", "ELM Distribution", "godis", size = "large", plotlyOutput("diselmchart")),
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
                          div(id="txtbox",textInput(inputId="text",label = "X-axis title", value = "Numbers of random DMIs")),
                          tags$style(type="text/css", "#txtbox {display: inline-block; max-width: 200px; }"),
                          div(id="txtbox", textInput("text2", label = "Y-axis title", value = "Frequency of random DMIs")),
                          div(id="txtbox",numericInput("xlimstart", label = "X-axis Start",0)),
                          div(id="txtbox",numericInput("xlimend", label = "X-axis End",500)),
                          tags$hr(),
                          tags$h4(tags$strong("Select Colors")),
                          
                          div(id="txtbox",colourInput("col", "Select bar colour", "deepskyblue1")),
                          div(id="txtbox",colourInput("col2", "Select background colour", "white")),
                          tags$hr(),
                          tags$h4(tags$strong("Select width/height to download plot as png")),
                          
                          div(id="txtbox",numericInput("width", label = "Width ", value = 1200)),
                          div(id="txtbox",numericInput("height", label = "Height ", value = 700))
                          
                          
                          
                          
                          
                      )),
                    tabPanel("Distribution of ELMs",
                             DT::dataTableOutput("diselmsdata"), tags$br(),tags$hr(),div(id="txtbox",actionButton("godis", "Interactive View"))
                    ),
                    tabPanel("Network",fluidPage(tags$br(), selectInput("selectlayout", label = "Select Layout",
                                                                        choices = list("Circle" = "layout_in_circle","Nice" = "layout_nicely", "Random" = "layout_randomly", "Piecewise" = "piecewise.layout", "Gem" = "layout.gem"),
                                                                        selected = "layout_in_circle"),
                                                 
                                                 
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
  #lowercase function
  lowername <- function(x) {
    colnames(x) <- tolower(colnames(x))
    x
    
  }

  #####################################################Domain-Motif Interactions####################################################
  #*********************************************************************************************************************************
  #Uploaded Data
  ####################################################
  
  inputDataMotif <-eventReactive(input$run, {
    MotifFile<-input$Motif
    validate(
      need(input$Motif != "", "Please select SLiM prediction file "), errorClass = "myClass"
    )
    
    #Read uploaded files
    Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
    
  })
  inputDataPPI <-eventReactive(input$run, {
    PPIFile<-input$PPI
    validate(
      need(input$PPI != "", "Please select PPI file "),errorClass = "myClass"
    )
    PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
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
  caption = tags$h4(tags$strong("SLiM Prediction File"))
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
    validate(
      need(input$Motif != "", "No SLiMs file was found Please upload file and 'Run' again "), errorClass = "myClass"
    )
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
    Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    
    
    #Motif-Domain Mapping
    #Rename the columns in two files
    names(Motif_NR) <- c("Seq", "Motif")
    names(Domain) <- c("Motif", "Domain")
    #Join/Merge two files based on Motif
    Join <- merge(Motif_NR, Domain, by="Motif")
    names(Join) <- c("Motif", "Seq", "Domain")  #Change header of the output file
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
    MotifFile<-input$Motif
    PPIFile<-input$PPI
    validate(
      need(input$Motif != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
    )
    validate(
      need(input$PPI != "", "PPI file is missing. Please upload and try again "), errorClass = "myClass"
    )
    validate(
      need(input$Motif != "", "Please select SLiMs file "), errorClass = "myClass"
    )
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
    Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    
    PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    #Rename the columns in two files
    names(Motif_NR) <- c("Seq", "Motif")
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
    #print(Uni_DMI)
    
    #PPI-DMI Mapping
    ########################################################################
    names(PPI2) <- c("mProtein", "dProtein")
    predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
    Uni_predDMIs <- unique(predDMI)
    names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
    predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
    print(predDMIs)
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
    MotifFile<-input$Motif
    PPIFile<-input$PPI
    validate(
      need(input$Motif != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
    )
    validate(
      need(input$PPI != "", "PPI file is missing. Please upload and try again "), errorClass = "myClass"
    )
    validate(
      need(input$Motif != "", "Please select SLiMs file "), errorClass = "myClass"
    )
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
    Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
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
      Datatype = factor(c("Motif","Domain","mProtien","dProtein")),
      Numbers = c(a,b,c,d)
    )
    
    p <- ggplot(data=statdmi, aes(x=Datatype, y=Numbers,fill=Datatype)) +
      geom_bar(colour="black", stat="identity") +
      guides(fill=FALSE)
    
    p <- ggplotly(p)
    #p
    
    #barplot(x, main="Statistics of DMIs", col = colors)
    
    #legend("topright",
    #legend = c(paste("Motif=",a),paste("Domain=",b),paste("Motif containing Proteins=",c),paste("Domain containing Proteins=",d)), fill = colors)
    
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
    MotifFile<-input$Motif
    PPIFile<-input$PPI
    validate(
      need(input$Motif != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
    )
    validate(
      need(input$PPI != "", "PPI file is missing. Please upload and try again "), errorClass = "myClass"
    )
    validate(
      need(input$Motif != "", "Please select SLiMs file "), errorClass = "myClass"
    )
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
    Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
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
    pvalueelm <- round(df_pred2$Frequency/nrow(df_pred),2)
    pvaluecol <- cbind(df_pred2,pvalueelm)
    names(pvaluecol) <- c("ELM", "Frequency", "Pvalue")
    GeneOntology <- merge(pvaluecol,GOterms, by="ELM")
    GeneOntology
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
    MotifFile<-input$Motif
    PPIFile<-input$PPI
    validate(
      need(input$Motif != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
    )
    validate(
      need(input$PPI != "", "PPI file is missing. Please upload and try again "), errorClass = "myClass"
    )
    validate(
      need(input$Motif != "", "Please select SLiMs file "), errorClass = "myClass"
    )
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
    Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
    Motif_NR<-unique(Motif)
    PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
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
    MotifFile<-input$Motif
    PPIFile<-input$PPI
    validate(
      need(input$Motif != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
    )
    validate(
      need(input$PPI != "", "PPI file is missing. Please upload and try again "), errorClass = "myClass"
    )
    validate(
      need(input$Motif != "", "Please select SLiMs file "), errorClass = "myClass"
    )
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
    #Read uploaded files
    Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
    Motif <- lowername(Motif)
    Motif <- Motif[, c("accnum","motif")]
    names(Motif) <- c("UniprotID","Motif")
     Motif_NR<-unique(Motif)
    PPI<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
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
    names(PPI) <- c("mProtein", "dProtein")
    names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
    predDMI <- merge(PPI, Uni_DMI, by= c("mProtein", "dProtein"))
    Uni_predDMIs <- unique(predDMI)
    names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
    predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
    #print(predDMIs)
    #Randomization/Permutations
    ##############################################################################
    
    PPIFile<-input$PPI

     PPI_data<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
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
    dirName <- paste0("RandomFiles_", strsplit(as.character(PPIFile$name), '.csv'))
    if(file.exists(paste0(dirName))){
      showNotification(paste0("RandomFiles folder already exists. Will be Loaded instead"),closeButton = TRUE, type = "warning")
      
      #checks if any of the 1000 random files is missing and creates it.
      for (l in 1:1000) {
        if(!file.exists(paste0(dirName,"/rPPI",l,".csv"))){
          
          permutation<-PermutationFunction(PPI_Matrix, k = length(PPI_Matrix))
          permutation2<-PermutationFunction(PPI_Matrix2, k = length(PPI_Matrix2))
          final_file<- c(paste(permutation,permutation2, sep = ":"))
          newCol1<-strsplit(as.character(final_file),':',fixed=TRUE)
          df<-data.frame(final_file,do.call(rbind, newCol1))
          subset = df[,c(2,3)]
          names(subset)<-c("mProtein", "dProtein")
          write.csv(subset, paste0(dirName,"/rPPI",l,".csv"), row.names = FALSE)
        }
      }
      #checks if randomNumber file exists and loads instead.
      if(file.exists(paste0(dirName,"/randomNumbers.csv"))){
        showNotification(paste0("randomNumbers file already exists. Loaded instead"),closeButton = TRUE, type = "warning")
        x<-read.csv(paste0(dirName,"/randomNumbers.csv"), sep = ",", header = FALSE)
      }
          
        
      
      else{
        for (i in 1:1000) {
          rPPI <- read.table(paste0(dirName,"/rPPI", i, ".csv"),
                             stringsAsFactors=FALSE, sep=",", strip.white=TRUE)
          names(rPPI)<-c("mProtein", "dProtein")
          names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
          DMI_rPPI <- merge(Uni_DMI, rPPI, by= c("mProtein", "dProtein"))
          Matches <- nrow(DMI_rPPI)
          print(Matches)
          output_randomNumbers <- write.table(Matches, paste0(dirName,"/randomNumbers.csv"), col.names = FALSE, append = TRUE, row.names = FALSE)

        }
      }
    }
    else{
      showNotification("Performing the randomizations", type = "message", duration = 5)
      dir.create(dirName)
    
        #for loop to create 1000 randomized files
        for (j in 1:1000) {
          permutation<-PermutationFunction(PPI_Matrix, k = length(PPI_Matrix))
          permutation2<-PermutationFunction(PPI_Matrix2, k = length(PPI_Matrix2))
          final_file<- c(paste(permutation,permutation2, sep = ":"))
          newCol1<-strsplit(as.character(final_file),':',fixed=TRUE)
          df<-data.frame(final_file,do.call(rbind, newCol1))
          subset = df[,c(2,3)]
          names(subset)<-c("mProtein", "dProtein")
          write.csv(subset, paste0(dirName,"/rPPI",j,".csv"), row.names = FALSE)
        
      }
      showNotification("1000 random PPI files have been created in RandomFiles folder", type = "message", duration = 5)
      showNotification("Now predicing DMIs from the random PPI data", type = "message", closeButton = TRUE,duration = 15)
      #rPPI-DMI Mapping
      #################################################################################
      for (i in 1:1000) {

        rPPI <- read.table(paste0(dirName,"/rPPI", i, ".csv"),
                           stringsAsFactors=FALSE, sep=",", strip.white=TRUE)
        names(rPPI)<-c("mProtein", "dProtein")
        names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
        DMI_rPPI <- merge(Uni_DMI, rPPI, by= c("mProtein", "dProtein"))
        Matches <- nrow(DMI_rPPI)
        print(Matches)
        
        output_randomNumbers <- write.table(Matches, paste0(dirName,"/randomNumbers.csv"), col.names = FALSE, append = TRUE, row.names = FALSE)

      }


    }

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
    MotifFile<-input$Motif
    PPIFile<-input$PPI
    validate(
      need(input$Motif != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
    )
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

    dirName <- paste0("RandomFiles_", strsplit(as.character(PPIFile$name), '.csv'))
    x<-read.csv(paste0(dirName,"/randomNumbers.csv"), sep = ",", header = FALSE)
    names(x) <- "values"
    bins<- seq(min(x), max(x), length.out = input$bins + 1)
    par(bg=input$col2)
    h<- hist(x$values, breaks=bins, col = input$col, border = 'black', main=input$text3, ylab=input$text2,
             xlab=input$text, xlim = c(0, max(x)+100), labels = input$barlabel, cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
    xfit <- seq(min(x$values),max(x$values),length=40)

    yfit <- dnorm(xfit, mean=mean(x$values),sd=sd(x$values))


    yfit<- yfit*diff(h$mids[1:2])*length(x$values)


    lines(xfit,yfit)
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
  output$summary <- renderUI({

    HTML(paste("<font color=\"#FF0000\"><b>Summary of Histogram</b></font>", pvalue, meanvalue, Escore, sep = '<hr/>'))

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
  #output$plot <- renderPlot({
   #  plotInput()
  #})

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
      MotifFile<-input$Motif
      PPIFile<-input$PPI
      validate(
        need(input$Motif != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
      )
      validate(
        need(input$PPI != "", "PPI file is missing. Please upload and try again "), errorClass = "myClass"
      )
      validate(
        need(input$Motif != "", "Please select SLiMs file "), errorClass = "myClass"
      )
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
      Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")
      Motif <- lowername(Motif)
      Motif <- Motif[, c("accnum","motif")]
      names(Motif) <- c("UniprotID","Motif")
      Motif_NR<-unique(Motif)
      PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
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
      png(file, width = input$width, height = input$height, units = "px", pointsize = 12)
      plotInput()
      dev.off()
    }
  )
})
}
# Run the application
shinyApp(ui = ui, server = server)


#*********************************************************************************************************
#*********************************************************************************************************
# Short Linear Motif Enrichment Analysis Server (SLiMEnrich)
# Developer: **Sobia Idrees**
# Version: 1.0.0
# Description: -SLiMEnrich predicts Domain Motif Interactions (DMIs) from Protein-Protein Interaction (PPI) data and analyzes enrichment through permutation test.
#              -This version doesn't save randomized files and matches file.
#              -Added better naming conventions (generic)
#*********************************************************************************************************
#*********************************************************************************************************
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
load_or_install(c("shiny", "ggplot2", "colourpicker", "shinyBS", "shinythemes", "DT", "shinyjs", "visNetwork", "igraph"))
#library(shiny)
#library(ggplot2)
#library(colourpicker)
#library(shinyBS)
#library(shinythemes)
#library(DT)
#library(shinyjs)
#library(visNetwork)
#library(igraph)
##############################
#GUI of the App
##############################
if (interactive()) {
  #navbar page with sidebar layout along with tabsets
  ui <- shinyUI(navbarPage(div(id="title", ("SLiMEnrich")), tabPanel("Domain-Motif Interactions", tags$head(
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
        bsModal("Hist", "Histogram", "go", size = "large", plotOutput("plot"),  downloadButton("downloadPlot", "Download")),
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
                      splitLayout(cellWidths = c("50%", "50%"), plotOutput("plotbar"))
                    )),
                    tabPanel("Histogram", fluidRow(
                      splitLayout(cellWidths = c("50%", "50%"), plotOutput("histogram"), htmlOutput("summary"))),
                      tags$hr(),
                      div(id="txtbox",actionButton("setting", "Settings")),
                      div(id="txtbox",actionButton("go", "Open plot in a separate window")),
                      
                      div(id="settings", sliderInput("bins",
                                                     "Number of bins",
                                                     min= 1,
                                                     max = 200,
                                                     value = 30),
                          tags$hr(),
                          checkboxInput("barlabel", label="Bar Labels", value = FALSE, width = NULL),
                          div(id="txtbox", textInput("text3", label = "Main title", value = "Distribution of random DMIs")),
                          div(id="txtbox",textInput(inputId="text",label = "X-axis title", value = "Numbers of random DMIs")),
                          tags$style(type="text/css", "#txtbox {display: inline-block; max-width: 200px; }"),
                          div(id="txtbox", textInput("text2", label = "Y-axis title", value = "Frequency of random DMIs")),
                          tags$hr(),
                          div(id="txtbox",colourInput("col", "Select bar colour", "deepskyblue1")),
                          div(id="txtbox",colourInput("col2", "Select background colour", "white"))
                          
                          
                          
                      )),
                    tabPanel("Network",fluidPage(tags$br(), selectInput("selectlayout", label = "Select Layout",
                                                                        choices = list("Nice" = "layout_nicely", "Circle" = "layout_in_circle", "Tree" = "layout_as_tree", "Davidson" = "layout.davidson.harel", "3D Grid" = "layout.grid.3d", "Sphere" = "layout_on_sphere", "Random" = "layout_randomly", "Piecewise" = "piecewise.layout", "Gem" = "layout.gem"),
                                                                        selected = "layout_nicely"), tags$img(src="key.png", align = "right"),
                                                 
                                                 
                                                 hr(),
                                                 
                                                 visNetworkOutput(outputId = "network",
                                                                  height = "1500px",
                                                                  width = "1500px")
                                                 
                                                 
                    )
                    )
                    
                    
        )
      )
    )),
    
    
    
    tabPanel("Getting Started"),tabPanel("Instructions", fluidPage(
      includeMarkdown("doc/instructions.Rmd")
    )), useShinyjs(),theme = shinytheme("cosmo"),
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
        dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")[,c(2:1)]
      }
      
      else{
        dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")[,c(2:1)]
      }
    })
    inputDataMotifDomain <-eventReactive(input$run, {
      MotifDomainFile<-input$MotifDomain
      if(is.null(MotifDomainFile)){
        Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")[,c(1:2)]
        
      }
      else{
        Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
      }
      
    })
    
    #shows the data table
    output$udata<-renderDataTable({
      
      inputDataMotif()
      
      
    })
    output$udata2<-renderDataTable({
      
      inputDataPPI()
      
    })
    output$udata3<-renderDataTable({
      DomainFile<-input$domain
      if(is.null(DomainFile)){
        return(NULL)
      }
      else{
        inputDatadomain()
      }
    })
    output$udata4<-renderDataTable({
      MotifDomainFile<-input$MotifDomain
      if(is.null(MotifDomainFile)){
        return(NULL)
      }
      else{
        inputDataMotifDomain()
      }
      
    })
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
        dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")[,c(2:1)]
      }
      
      else{
        dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")[,c(2:1)]
      }
      MotifDomainFile<-input$MotifDomain
      if(is.null(MotifDomainFile)){
        Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")[,c(1:2)]
        
      }
      else{
        Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
      }
      #Read uploaded files
      Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")[,c('AccNum','Motif')]
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
        dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")[,c(2:1)]
      }
      
      else{
        dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")[,c(2:1)]
      }
      MotifDomainFile<-input$MotifDomain
      if(is.null(MotifDomainFile)){
        Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")[,c(1:2)]
        
      }
      else{
        Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
      }
      #Read uploaded files
      Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")[,c('AccNum','Motif')]
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
      
      #vhPPI-DMI Mapping
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
        dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")[,c(2:1)]
      }
      
      else{
        dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")[,c(2:1)]
      }
      MotifDomainFile<-input$MotifDomain
      if(is.null(MotifDomainFile)){
        Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")[,c(1:2)]
        
      }
      else{
        Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
      }
      
      #Read uploaded files
      Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")[,c('AccNum','Motif')]
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
      #vhPPI-DMI Mapping
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
      #Label names for the chart
      labels <- c("Motifs", "Domains", "Motif containing Proteins", "Domain containing Proteins")
      #Created bar char of the unique values
      barplot(x, main="Statistics of DMIs", col = colors)
      legend("topright",
             legend = c(paste("Motif=",a),paste("Domain=",b),paste("Motif containing Proteins=",c),paste("Domain containing Proteins=",d)), fill = colors)
      
    })
    output$plotbar <- renderPlot({
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

  #####################################################Enrichment Analysis####################################################
  #*********************************************************************************************************************************


  #*************************************************************************************

  #Step 4: rPPI-DMI Mapping
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
        dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")[,c(2:1)]
      }
      
      else{
        dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")[,c(2:1)]
      }
      MotifDomainFile<-input$MotifDomain
      if(is.null(MotifDomainFile)){
        Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")[,c(1:2)]
        
      }
      else{
        Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
      }
      observe({
        show(id = "go2")
      })
      #Read uploaded files
      Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")[,c('AccNum','Motif')]
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
    #dir.create("RandomFiles")
    
    #dirName <- paste0("RandomFiles_", strsplit(as.character(PPIFile$name), '.csv'))
    #dir.create(dirName)
    #for loop to create 1000 randomized files
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
    
    #output_randomNumbers <- write.csv(m, paste0("randomNumbers.csv"), row.names = FALSE)
    m
  })
  #*************************************************************************************

  #Step 5: Histogram
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
          progress$set(value = value, detail = "This may take a while!!")
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
        dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")[,c(2:1)]
      }
      
      else{
        dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")[,c(2:1)]
      }
      MotifDomainFile<-input$MotifDomain
      if(is.null(MotifDomainFile)){
        Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")[,c(1:2)]
        
      }
      else{
        Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
      }
      
    #dirName <- paste0("RandomFiles_", strsplit(as.character(PPIFile$name), '.csv'))
    x <- rPPIDMI()
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
    FDR <- paste0("<b>False Discrovery Rate is: </b>", round(mean(x$values)/nrow(predictedDMIs()),2))
    output$summary <- renderUI({
      
      HTML(paste("<font color=\"#FF0000\"><b>Summary of Histogram</b></font>", pvalue, meanvalue, FDR, sep = '<hr/>'))
      
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

  #Step 6: DMI Network
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
        dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")[,c(2:1)]
      }
      
      else{
        dProtein<-read.csv(DomainFile$datapath,header=TRUE,sep=",")[,c(2:1)]
      }
      MotifDomainFile<-input$MotifDomain
      if(is.null(MotifDomainFile)){
        Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")[,c(1:2)]
        
      }
      else{
        Domain<-read.csv(MotifDomainFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
      }
      #Read uploaded files
      Motif<-read.csv(MotifFile$datapath,header=TRUE,sep=",")[,c('AccNum','Motif')]
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
      
      #vhPPI-DMI Mapping
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
        visExport(g6, type = "png", name = "export-network",float = "left", label = "Save network", background = "white", style= "")
        
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
      png(file, width = 1200, height = 700, units = "px", pointsize = 12)
      plotInput()
      dev.off()
    }
  )
  })
}
# Run the application
shinyApp(ui = ui, server = server)


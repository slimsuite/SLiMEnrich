#*********************************************************************************************************
#*********************************************************************************************************
# Short Linear Motif Enrichment Analysis App (SLiMEnrich)
# Developer: **Sobia Idrees**
# Version: 1.0.1
# Description: SLiMEnrich predicts Domain Motif Interactions (DMIs) from Protein-Protein Interaction (PPI) data and analyzes enrichment through permutation test.
#*********************************************************************************************************
#*********************************************************************************************************
##############################
#Version History
##############################
#V1.0.1 - Added code for checking whether packages installed. (Removes manual step)
##############################


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

                                                                       fileInput("ELM","Select SLiMProb.occ file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
                                                                       actionButton("run", "Run", width = "100px"),
                                                                       div(id="fileuploads",checkboxInput("uploadelms",label = "Upload Domain and/or ELM.Pfam Files", value = FALSE)),
                                                                       div(id="uploadelm",  fileInput("domain","Select Domain file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
                                                                           fileInput("ELMPfam","Select ELM-Pfam file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')))
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
  observeEvent(input$uploadelms, {
    toggle(id = "uploadelm", anim = TRUE)
  })

  #####################################################Domain-Motif Interactions####################################################
  #*********************************************************************************************************************************
  #Uploaded Data
  ####################################################

  inputDataELM <-eventReactive(input$run, {
    ELMFile<-input$ELM
    validate(
      need(input$ELM != "", "Please select SLiMs file "), errorClass = "myClass"
    )

    #Read uploaded files
    ELM<-read.csv(ELMFile$datapath,header=TRUE,sep=",")

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
    PfamFile<-input$domain
    if(is.null(PfamFile)){
      hProtein<-read.csv("data/Pfams.csv",header=TRUE,sep=",")[,c('pfam','accnum')]
    }

    else{
      hProtein<-read.csv(PfamFile$datapath,header=TRUE,sep=",")[,c('pfam','accnum')]
    }
  })
  inputDataELMPfam <-eventReactive(input$run, {
    ELMPfamFile<-input$ELMPfam
    if(is.null(ELMPfamFile)){
      Pfam<-read.csv("data/ELM.Pfam.tsv",header=TRUE,sep="\t")[,c(1:2)]

    }
    else{
      Pfam<-read.csv(ELMPfamFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
    }

  })

  #shows the data table
  output$udata<-renderDataTable({

    inputDataELM()


  })
  output$udata2<-renderDataTable({

    inputDataPPI()

  })
  output$udata3<-renderDataTable({
    PfamFile<-input$domain
    if(is.null(PfamFile)){
      return(NULL)
    }
    else{
    inputDatadomain()
    }
  })
  output$udata4<-renderDataTable({
    ELMPfamFile<-input$ELMPfam
    if(is.null(ELMPfamFile)){
      return(NULL)
    }
    else{
    inputDataELMPfam()
    }

  })
  #*************************************************************************************
  #Step 1: Potential DMIs
  ####################################################
  potentialDMIs <-eventReactive(input$run, {
    #File upload check
    ELMFile<-input$ELM
    validate(
      need(input$ELM != "", "No SLiMs file was found Please upload file and 'Run' again "), errorClass = "myClass"
    )
    #File upload check
    PfamFile<-input$domain
    if(is.null(PfamFile)){
      hProtein<-read.csv("data/Pfams.csv",header=TRUE,sep=",")[,c('pfam','accnum')]
    }

    else{
      hProtein<-read.csv(PfamFile$datapath,header=TRUE,sep=",")[,c('pfam','accnum')]
    }
    ELMPfamFile<-input$ELMPfam
    if(is.null(ELMPfamFile)){
      Pfam<-read.csv("data/ELM.Pfam.tsv",header=TRUE,sep="\t")[,c(1:2)]

  }
    else{
      Pfam<-read.csv(ELMPfamFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
  }
    #Read uploaded files
    ELM<-read.csv(ELMFile$datapath,header=TRUE,sep=",")[,c('AccNum','Motif')]
    ELM_NR<-unique(ELM)


#ELM-Pfam Mapping
    #Rename the columns in two files
    names(ELM_NR) <- c("Seq", "ELM")
    names(Pfam) <- c("ELM", "Pfam")
    #Join/Merge two files based on ELM
    Join <- merge(ELM_NR, Pfam, by="ELM")
    #print(Join)
    names(Join) <- c("ELM", "Seq", "Pfams")  #Change header of the output file
#Pfam-hProtein Mapping
    #Load results from the previous code)
    names(hProtein) <- c("Pfams", "hProteins")
    #joined both files based on Pfams
    DMI <- merge(Join, hProtein,by="Pfams")
    #Filtered unique DMIs
    Uni_DMI <- unique(DMI)
    #Named the header of output file
    names(Uni_DMI) <- c("Pfam", "ELM", "Protein1", "Protein2")
    Uni_DMI <- Uni_DMI[, c("Protein1","ELM", "Pfam", "Protein2")]
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
    ELMFile<-input$ELM
    PPIFile<-input$PPI
    validate(
      need(input$ELM != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
    )
    validate(
      need(input$PPI != "", "PPI file is missing. Please upload and try again "), errorClass = "myClass"
    )
    validate(
      need(input$ELM != "", "Please select SLiMs file "), errorClass = "myClass"
    )
    #File upload check
    PfamFile<-input$domain
    if(is.null(PfamFile)){
      hProtein<-read.csv("data/Pfams.csv",header=TRUE,sep=",")[,c('pfam','accnum')]
    }

    else{
      hProtein<-read.csv(PfamFile$datapath,header=TRUE,sep=",")[,c('pfam','accnum')]
    }
    ELMPfamFile<-input$ELMPfam
    if(is.null(ELMPfamFile)){
      Pfam<-read.csv("data/ELM.Pfam.tsv",header=TRUE,sep="\t")[,c(1:2)]

    }
    else{
      Pfam<-read.csv(ELMPfamFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
    }
    #Read uploaded files
    ELM<-read.csv(ELMFile$datapath,header=TRUE,sep=",")[,c('AccNum','Motif')]
    ELM_NR<-unique(ELM)

    PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    #Rename the columns in two files
    names(ELM_NR) <- c("Seq", "ELM")
    names(Pfam) <- c("ELM", "Pfam")
    #Join/Merge two files based on ELM
    Join <- merge(ELM_NR, Pfam, by="ELM")
    #print(Join)
    names(Join) <- c("ELM", "Seq", "Pfams")  #Change header of the output file
    #Load vORF_ELM_Pfam file (result file from the previous code)
    names(hProtein) <- c("Pfams", "hProteins")
    #joined both files based on Pfams
    DMI <- merge(Join, hProtein,by="Pfams")
    #Filtered unique DMIs
    Uni_DMI <- unique(DMI)
    #Named the header of output file
    names(Uni_DMI) <- c("Pfam", "ELM", "vProtein", "hProtein")
    #print(Uni_DMI)

    #vhPPI-DMI Mapping
    ########################################################################
    names(PPI2) <- c("vProtein", "hProtein")
    predDMI <- merge(PPI2, Uni_DMI, by= c("vProtein", "hProtein"))
    Uni_predDMIs <- unique(predDMI)
    names(Uni_predDMIs) <- c("Protein1", "Protein2", "Pfam", "ELM")
    predDMIs <- Uni_predDMIs[, c("Protein1","ELM", "Pfam", "Protein2")]
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
    ELMFile<-input$ELM
    PPIFile<-input$PPI
    validate(
      need(input$ELM != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
    )
    validate(
      need(input$PPI != "", "PPI file is missing. Please upload and try again "), errorClass = "myClass"
    )
    validate(
      need(input$ELM != "", "Please select SLiMs file "), errorClass = "myClass"
    )
    #File upload check
    PfamFile<-input$domain
    if(is.null(PfamFile)){
      hProtein<-read.csv("data/Pfams.csv",header=TRUE,sep=",")[,c('pfam','accnum')]
    }

    else{
      hProtein<-read.csv(PfamFile$datapath,header=TRUE,sep=",")[,c('pfam','accnum')]
    }
    ELMPfamFile<-input$ELMPfam
    if(is.null(ELMPfamFile)){
      Pfam<-read.csv("data/ELM.Pfam.tsv",header=TRUE,sep="\t")[,c(1:2)]

    }
    else{
      Pfam<-read.csv(ELMPfamFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
    }

      #Read uploaded files
      ELM<-read.csv(ELMFile$datapath,header=TRUE,sep=",")[,c('AccNum','Motif')]
      ELM_NR<-unique(ELM)
      PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
      #Rename the columns in two files
      names(ELM_NR) <- c("Seq", "ELM")
      names(Pfam) <- c("ELM", "Pfam")
      #Join/Merge two files based on ELM
      Join <- merge(ELM_NR, Pfam, by="ELM")
      #print(Join)
      names(Join) <- c("ELM", "Seq", "Pfams")  #Change header of the output file
      #Load vORF_ELM_Pfam file (result file from the previous code)
      names(hProtein) <- c("Pfams", "hProteins")
      #joined both files based on Pfams
      DMI <- merge(Join, hProtein,by="Pfams")
      #Filtered unique DMIs
      Uni_DMI <- unique(DMI)
      #Named the header of output file
      names(Uni_DMI) <- c("Pfam", "ELM", "vProtein", "hProtein")
      #print(Uni_DMI)
      #vhPPI-DMI Mapping
      ########################################################################
      names(PPI2) <- c("vProtein", "hProtein")
      predDMI <- merge(PPI2, Uni_DMI, by= c("vProtein", "hProtein"))
      Uni_predDMIs <- unique(predDMI)
      names(Uni_predDMIs) <- c("Protein1", "Protein2", "Pfam", "ELM")
      predDMIs <- Uni_predDMIs[, c("Protein1","ELM", "Pfam", "Protein2")]
      #print(predDMIs)
      colors=c("cadetblue1", "deepskyblue2", "blue", "darkblue")
      #Select unique ELM
      uniq_elm <- unique(predDMI$ELM)
      a <- length(uniq_elm)
      #Select unique Pfam
      uniq_pfam <- unique(predDMI$Pfam)
      b <- length(uniq_pfam)
      #Select unique vORF
      uniq_vorf <- unique(predDMI$vProtein)
      c <- length(uniq_vorf)
      #Select unique hproteins
      uniq_hprotein <- unique(predDMI$hProtein)
      d <- length(uniq_hprotein)
      uniq_count<-c(a = length(uniq_elm), b = length(uniq_pfam), c = length(uniq_vorf), d = length(uniq_hprotein))
      #Create pie chart
      x <- c(a,b,c,d)
      #Label names for the chart
      labels <- c("ELMs", "Pfams", "ELM containing Proteins", "Domain containing Proteins")
      #Created bar char of the unique values
      barplot(x, main="Statistics of DMIs", col = colors)
      legend("topright",
             legend = c(paste("ELM=",a),paste("Pfam=",b),paste("ELM containing Proteins=",c),paste("Domain containing Proteins=",d)), fill = colors)

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
    ELMFile<-input$ELM
    PPIFile<-input$PPI
    validate(
      need(input$ELM != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
    )
    validate(
      need(input$PPI != "", "PPI file is missing. Please upload and try again "), errorClass = "myClass"
    )
    validate(
      need(input$ELM != "", "Please select SLiMs file "), errorClass = "myClass"
    )
    #File upload check
    PfamFile<-input$domain
    if(is.null(PfamFile)){
      hProtein<-read.csv("data/Pfams.csv",header=TRUE,sep=",")[,c('pfam','accnum')]
    }

    else{
      hProtein<-read.csv(PfamFile$datapath,header=TRUE,sep=",")[,c('pfam','accnum')]
    }
    ELMPfamFile<-input$ELMPfam
    if(is.null(ELMPfamFile)){
      Pfam<-read.csv("data/ELM.Pfam.tsv",header=TRUE,sep="\t")[,c(1:2)]

    }
    else{
      Pfam<-read.csv(ELMPfamFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
    }
    observe({
      show(id = "go2")
    })
    #Read uploaded files
    ELM<-read.csv(ELMFile$datapath,header=TRUE,sep=",")[,c('AccNum','Motif')]
    ELM_NR<-unique(ELM)
    PPI<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    #Rename the columns in two files
    names(ELM_NR) <- c("Seq", "ELM")
    names(Pfam) <- c("ELM", "Pfam")
    #Join/Merge two files based on ELM
    Join <- merge(ELM_NR, Pfam, by="ELM")
    #print(Join)
    names(Join) <- c("ELM", "Seq", "Pfams")  #Change header of the output file
    #Load vORF_ELM_Pfam file (result file from the previous code)
    names(hProtein) <- c("Pfams", "hProteins")
    #joined both files based on Pfams
    DMI <- merge(Join, hProtein,by="Pfams")
    #Filtered unique DMIs
    Uni_DMI <- unique(DMI)
    #Named the header of output file
    names(PPI) <- c("vProtein", "hProtein")
    names(Uni_DMI) <- c("Pfam", "ELM", "vProtein", "hProtein")
    predDMI <- merge(PPI, Uni_DMI, by= c("vProtein", "hProtein"))
    Uni_predDMIs <- unique(predDMI)
    names(Uni_predDMIs) <- c("Protein1", "Protein2", "Pfam", "ELM")
    predDMIs <- Uni_predDMIs[, c("Protein1","ELM", "Pfam", "Protein2")]
    #print(predDMIs)
    #Randomization/Permutations
    ##############################################################################
    PPIFile<-input$PPI

     PPI_data<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
    names(PPI_data) <- c("vORF", "hProtein")
    PPI_Matrix<-matrix(data = PPI_data$vORF)
    PPI_Matrix2<-matrix(data = PPI_data$hProtein)
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
      showNotification(paste0("RandomFiles folder already exists. Loaded instead"),closeButton = TRUE, type = "warning")
      if(file.exists(paste0(dirName,"/randomNumbers.csv"))){
        showNotification(paste0("randomNumbers file already exists. Loaded instead"),closeButton = TRUE, type = "warning")
        x<-read.csv(paste0(dirName,"/randomNumbers.csv"), sep = ",", header = FALSE)
      }
      else{
        for (i in 1:1000) {
          rPPI <- read.table(paste0(dirName,"/rPPI", i, ".csv"),
                             stringsAsFactors=FALSE, sep=",", strip.white=TRUE)
          names(rPPI)<-c("vProtein", "hProtein")
          names(Uni_DMI) <- c("Pfam", "ELM", "vProtein", "hProtein")
          DMI_rPPI <- merge(Uni_DMI, rPPI, by= c("vProtein", "hProtein"))
          Matches <- nrow(DMI_rPPI)
          print(Matches)
          output_randomNumbers <- write.table(Matches, paste0(dirName,"/randomNumbers.csv"), col.names = FALSE, append = TRUE, row.names = FALSE)

        }
      }
    }
    else{
      dir.create(dirName)
        #for loop to create 1000 randomized files
        for (j in 1:1000) {
          permutation<-PermutationFunction(PPI_Matrix, k = length(PPI_Matrix))
          permutation2<-PermutationFunction(PPI_Matrix2, k = length(PPI_Matrix2))
          final_file<- c(paste(permutation,permutation2, sep = ":"))
          newCol1<-strsplit(as.character(final_file),':',fixed=TRUE)
          df<-data.frame(final_file,do.call(rbind, newCol1))
          subset = df[,c(2,3)]
          names(subset)<-c("vProtein", "hProtein")
          write.csv(subset, paste0(dirName,"/rPPI",j,".csv"), row.names = FALSE)

        }
      showNotification("1000 random PPI files have been created in RandomFiles folder", type = "message", duration = 5)
      showNotification("Now predicing DMIs from the random PPI data", type = "message", closeButton = TRUE,duration = 10)
      #rPPI-DMI Mapping
      #################################################################################
      for (i in 1:1000) {

        rPPI <- read.table(paste0(dirName,"/rPPI", i, ".csv"),
                           stringsAsFactors=FALSE, sep=",", strip.white=TRUE)
        names(rPPI)<-c("vProtein", "hProtein")
        names(Uni_DMI) <- c("Pfam", "ELM", "vProtein", "hProtein")
        DMI_rPPI <- merge(Uni_DMI, rPPI, by= c("vProtein", "hProtein"))
        Matches <- nrow(DMI_rPPI)
        print(Matches)
        output_randomNumbers <- write.table(Matches, paste0(dirName,"/randomNumbers.csv"), col.names = FALSE, append = TRUE, row.names = FALSE)

      }


    }

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
    ELMFile<-input$ELM
    PPIFile<-input$PPI
    validate(
      need(input$ELM != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
    )
    #File upload check
    PfamFile<-input$domain
    if(is.null(PfamFile)){
      hProtein<-read.csv("data/Pfams.csv",header=TRUE,sep=",")[,c('pfam','accnum')]
    }

    else{
      hProtein<-read.csv(PfamFile$datapath,header=TRUE,sep=",")[,c('pfam','accnum')]
    }
    ELMPfamFile<-input$ELMPfam
    if(is.null(ELMPfamFile)){
      Pfam<-read.csv("data/ELM.Pfam.tsv",header=TRUE,sep="\t")[,c(1:2)]

    }
    else{
      Pfam<-read.csv(ELMPfamFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
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
      ELMFile<-input$ELM
      PPIFile<-input$PPI
      validate(
        need(input$ELM != "" || input$PPI != "", "Either SLiMs or PPI file is missing. Please upload files and try again "), errorClass = "myClass"
      )
      validate(
        need(input$PPI != "", "PPI file is missing. Please upload and try again "), errorClass = "myClass"
      )
      validate(
        need(input$ELM != "", "Please select SLiMs file "), errorClass = "myClass"
      )
      #File upload check
      PfamFile<-input$domain
      if(is.null(PfamFile)){
        hProtein<-read.csv("data/Pfams.csv",header=TRUE,sep=",")[,c('pfam','accnum')]
      }

      else{
        hProtein<-read.csv(PfamFile$datapath,header=TRUE,sep=",")[,c('pfam','accnum')]
      }
      ELMPfamFile<-input$ELMPfam
      if(is.null(ELMPfamFile)){
        Pfam<-read.csv("data/ELM.Pfam.tsv",header=TRUE,sep="\t")[,c(1:2)]

      }
      else{
        Pfam<-read.csv(ELMPfamFile$datapath,header=TRUE,sep="\t")[,c(1:2)]
      }
      #Read uploaded files
      ELM<-read.csv(ELMFile$datapath,header=TRUE,sep=",")[,c('AccNum','Motif')]
      ELM_NR<-unique(ELM)
      PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
      #Rename the columns in two files
      names(ELM_NR) <- c("Seq", "ELM")
      names(Pfam) <- c("ELM", "Pfam")
      #Join/Merge two files based on ELM
      Join <- merge(ELM_NR, Pfam, by="ELM")
      #print(Join)
      names(Join) <- c("ELM", "Seq", "Pfams")  #Change header of the output file
      #Load vORF_ELM_Pfam file (result file from the previous code)
      names(hProtein) <- c("Pfams", "hProteins")
      #joined both files based on Pfams
      DMI <- merge(Join, hProtein,by="Pfams")
      #Filtered unique DMIs
      Uni_DMI <- unique(DMI)
      #Named the header of output file
      names(Uni_DMI) <- c("Pfam", "ELM", "vProtein", "hProtein")
      #print(Uni_DMI)

      #vhPPI-DMI Mapping
      ########################################################################
      names(PPI2) <- c("vProtein", "hProtein")
      predDMI <- merge(PPI2, Uni_DMI, by= c("vProtein", "hProtein"))

      Uni_predDMIs <- unique(predDMI)
      names(Uni_predDMIs) <- c("Protein1", "Protein2", "Pfam", "ELM")
      predDMIs <- Uni_predDMIs[, c("Protein1","ELM", "Pfam", "Protein2")]
      if(nrow(predDMIs) != 0){
      #Network

      first <- predDMI[,c("vProtein","ELM")]
      g <- graph.data.frame(first, directed = F)
      #igraph.options(plot.layout=layout.graphopt, vertex.size=10)
      V(g)$color <- ifelse(V(g)$name %in% predDMI[,1], "#A93226", "#F7DC6F")
      V(g)$shape <- ifelse(V(g)$name %in% predDMI[,1], "square", "box")
      #plot(g,  edge.color="orange")
      #visIgraph(g)
      ##print(first)
      #V(g)$color <- "red"
      #Second
      second <- predDMI[,c("ELM","Pfam")]
      g2 <- graph.data.frame(second, directed = F)
      #igraph.options(plot.layout=layout.graphopt, vertex.size=10)
      V(g2)$color <- ifelse(V(g2)$name %in% predDMI[,1], "#9B59B6", "#D35400")
      V(g2)$shape <- ifelse(V(g2)$name %in% predDMI[,1], "box", "circle")
      #plot(g2,  edge.color="orange")
      #print(second)
      #visIgraph(g2)
      #V(g2)$color <- "green"
      #Third

      third <- predDMI[,c("Pfam","hProtein")]
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


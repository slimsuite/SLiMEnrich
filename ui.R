#*********************************************************************************************************
#*********************************************************************************************************
# Please see main.R for App version, history and license information.
source("main.R")
#*********************************************************************************************************
#*********************************************************************************************************
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
                  #info{
                  background-color: white;
                  color: black;
                  font-weight: bold;
                  font-size: 10;
                  }
                  
                  " ))
  ),
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      tags$div(class="note", checked=NA,
               tags$h4(paste(info$apptitle,"Version",info$version))#,
               #tags$b("Note: To analyse example dataset, press 'Load Data' without uploading any files."),
      ),
      
      fileInput("PPI","Select Interaction file:",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
      # PPI fields
      tags$div(class="ppifields", checked=NA,
               textInput(inputId="ppimprotein",label = "Motif-containing protein column", value = "mProtein"),
               textInput(inputId="ppidprotein",label = "Domain-containing protein column", value = "dProtein")
      ),
      
      
      actionButton("run", "Load data", width = "100px"),
      hr(),
      
      prettyRadioButtons(inputId = "DMIStrategy",
                         label = "DMI Strategy", icon = icon("check"),
                         choices = c("Link mProteins directly to dProteins (ELMi-Protein)" = "elmiprot", 
                                     "Link Motif classes directly to dProteins (ELMc-Protein)" = "elmcprot",
                                     "Link Motif classes to binding domains (ELMc-Domain)" = "elmcdom"),
                         selected = "elmcprot",
                         animation = "pulse", status = "warning"),
      #hr(),
      div(id="hidehelpmd",prettyCheckbox("hidehelp",label = tags$b("Hide tab info text"), value = FALSE, status = "info",
                                          icon = icon("check"),
                                          animation = "pulse")),
      div(id="fileuploads",prettyCheckbox("uploadmotifs",label = tags$b("Upload Additional Files"), value = FALSE, status = "info",
                                          icon = icon("check"),
                                          animation = "pulse")),
      #hr(),
      
      div(id="uploadmotif", 
          fileInput("MotifDomain","Select Motif-Domain file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
          conditionalPanel(
            condition = "input.SLiMrunid == false",
            fileInput("Motif","Select Motif file (mProtein-Motif, e.g. ELM or SLiMProb)",accept=c('text/csv','text/comma-separated-values,text/plain','csv'))
          ),
          div(id = "slimrun", textInput("SLiMRun", label = "SLiMProb JobID:", value = "")),
          prettyCheckbox("SLiMrunid", label = "Provide SLiMProb Job ID (Replaces Motif File)", status = "default",
                         icon = icon("check"),
                         animation = "pulse"),
          fileInput("domain","Select Domain file (Domain-dProtein)",accept=c('text/csv','text/comma-separated-values,text/plain','csv'))
      ),
      #div (id = "note", "Note: To analyse example dataset, press 'Load Data' without uploading any files"),
      #hr(),
      div(id="advsettings", 
          # DMI fields
          tags$div(class="dmifields", checked=NA,
                   #tags$hr(),
                   tags$h4("DMI File:"),
                   tags$p("Motifs and interacting Domain fields"),
                   # Default fields are from the ELM interactions table
                   textInput(inputId="dmimotif",label = "Motif column", value = "Elm"),
                   textInput(inputId="dmidomain",label = "Domain column", value = "interactorDomain")
          ),
          # Motif fields
          tags$div(class="motfields", checked=NA,
                   #tags$hr(),
                   tags$h4("Motifs File:"),
                   tags$p("Motif-containing proteins and their motifs"),
                   # Default fields are from the ELM instances table, reformatted to match SLiMProb
                   textInput(inputId="motifmprotein",label = "Motif file mProtein column", value = "AccNum"),
                   textInput(inputId="motifmotif",label = "Motif file Motif column", value = "Motif")
          ),
          # Domain fields
          tags$div(class="domfields", checked=NA,
                   #tags$hr(),
                   tags$h4("Domains File:"),
                   tags$p("Domain-containing proteins and their domains"),
                   # Default fields are from the Uniprot Pfam domain table
                   textInput(inputId="domaindomain",label = "Domain file domain column", value = "pfam"),
                   textInput(inputId="domaindprotein",label = "Domain file dProtein column", value = "accnum")
          )
      ),
      # Randomisation fields
      tags$div(class="header", checked=NA,
               #tags$hr(),
               tags$h4("Randomisation settings:"),
               numericInput("shufflenum", label = "Number of randomisations",1000,step=100,min=100)
      ),
      div (id = "update", paste0("Last updated:", info$lastedit))
    ),
    
    # MainPanel
    mainPanel(
      #Creates a seperate window (pop up window)
      bsModal("DisE", "ELM Distribution", "godis", size = "large", plotlyOutput("diselmchart")),
      bsModal("DidsE", "Domain Distribution", "godisd", size = "large", plotlyOutput("disdomchart")),
      #Tab view
      tabsetPanel(type="tabs",
                  tabPanel("Uploaded Data",
                           #htmlOutput("docs_tables"),
                           conditionalPanel(
                             condition = "input.hidehelp == false",
                             includeMarkdown("doc/tabs/uploaded.md")
                           ),
                           div(id="fullfilecheck",prettyCheckbox("parseddata",label = tags$b("Show parsed data columns"), value = FALSE, status = "info",
                                                                 icon = icon("check"),
                                                                 animation = "pulse")),
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%", "50%", "50%"), DT::dataTableOutput("udata2"), DT::dataTableOutput("udata")), DT::dataTableOutput("udata4"), DT::dataTableOutput("udata3")
                           )
                  ),
                  
                  tabPanel("Potential DMIs",
                           div(id="nrpotdmicheck",prettyCheckbox("nrpotdmi",label = tags$b("Show NR potential DMI"), value = FALSE, status = "info",
                                                                 icon = icon("check"),
                                                                 animation = "pulse")),
                           DT::dataTableOutput("data"),
                           tags$hr(),
                           downloadButton('downloadDMI', 'Download')
                  ),
                  
                  tabPanel("Predicted DMIs", 
                           div(id="nrdmicheck",prettyCheckbox("nrdmi",label = tags$b("Show NR predicted DMI"), value = FALSE, status = "info",
                                                              icon = icon("check"),
                                                              animation = "pulse")),
                           DT::dataTableOutput("PredDMIs"),tags$hr(),downloadButton('downloadpredDMI', 'Download')
                  ),
                  
                  tabPanel("Statistics", fluidRow(
                    splitLayout(cellWidths = c("75%", "25%"), plotlyOutput("plotbar"))
                  )),
                  tabPanel("Histogram", fluidRow(
                    splitLayout(cellWidths = c("50%", "50%"), plotOutput("histogram"), htmlOutput("summary"))),
                    tags$hr(),
                    div(id="txtbox",actionButton("setting", "Settings")),
                    div(id="txtbox",downloadButton("downloadPlot", "Download")),
                    
                    div(id="settings", 
                        #sliderInput("bins", "Number of bins", min= 1, max = 200, value = 30),
                        sliderInput("binwidth", "Width of bins", min= 1, max = 100, value = 1),
                        tags$hr(),
                        tags$h4(tags$strong("Select labels")),
                        
                        checkboxInput("barlabel", label="Bar Labels", value = FALSE, width = NULL),
                        div(id="txtbox", textInput("text3", label = "Main title", value = "Distribution of random DMIs")),
                        div(id="txtbox",textInput(inputId="text",label = "X-axis title", value = "Number of random DMIs")),
                        tags$style(type="text/css", "#txtbox {display: inline-block; max-width: 200px; }"),
                        div(id="txtbox", textInput("text2", label = "Y-axis title", value = "Frequency of random DMIs")),
                        div(id="txtbox",numericInput("xlimstart", label = "X-axis Start",0)),
                        div(id="txtbox",numericInput("xlimend", label = "Extend X-axis End",0)),
                        tags$hr(),
                        tags$h4(tags$strong("Select Colors")),
                        
                        div(id="txtbox",colourpicker::colourInput("col", "Select bar colour", "firebrick3")),
                        div(id="txtbox",colourpicker::colourInput("col2", "Select background colour", "white")),
                        tags$hr(),
                        tags$h4(tags$strong("Select width/height to download plot as png")),
                        
                        div(id="txtbox",numericInput("width", label = "Width ", value = 2600)),
                        div(id="txtbox",numericInput("height", label = "Height ", value = 1600))
                        
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
                  ),
                  tabPanel("Help",
                           includeMarkdown("README.md")     
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
# End of GUI code.
##############################
#*********************************************************************************************************
#*********************************************************************************************************


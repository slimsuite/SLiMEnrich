#*********************************************************************************************************
#*********************************************************************************************************
# Short Linear Motif Enrichment Analysis App (SLiMEnrich)
# Developer: **Sobia Idrees**
# Version: 1.0.8
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
##############################
#SLiMEnrich is free software: you can redistribute it and/or modify
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
package_names = c("shiny", "ggplot2", "colourpicker", "shinyBS", "shinythemes", "DT", "shinyjs", "visNetwork", "igraph","markdown","plotly", "plyr")
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
      
      fileInput("Motif","Select SLiM prediction file (e.g. SLiMProb)",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
      actionButton("run", "Run", width = "100px"),
      div(id="fileuploads",checkboxInput("uploadmotifs",label = "Upload Domain and/or Motif-Domain Files", value = FALSE)),
      div(id="uploadmotif",  fileInput("domain","Select Domain file",accept=c('text/csv','text/comma-separated-values,text/plain','csv')),
          fileInput("MotifDomain","Select Motif-Domain file",accept=c('text/csv','text/comma-separated-values,text/plain','csv'))),
      div (id = "note", "Note: To analyze example dataset, press 'Run' without uploading any files"),
      hr(),
      div (id = "update", "Last updated: 06-Sep-2017")
    ),
    
    # MainPanel
    mainPanel(
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



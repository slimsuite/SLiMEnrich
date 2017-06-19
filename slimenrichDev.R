#*********************************************************************************************************
#*********************************************************************************************************
# Short Linear Motif Enrichment Analysis App (SLiMEnrich)
# Developer: **Sobia Idrees**
# Version: 1.0.5
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
#V1.0.5 - Added a new tab to show distribution of ELMs in the predicted DMI dataset in tabular as well as interactive view.
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
load_or_install(c("shiny", "ggplot2", "colourpicker", "shinyBS", "shinythemes", "DT", "shinyjs", "visNetwork", "igraph","markdown","plotly", "plyr"))

##############################
#GUI of the App
##############################
if (interactive()) {
sourceDir <- getSrcDirectory(function(dummy) {dummy})
source(paste(sourceDir, "/ui.R", sep=""))

##############################
# Server logic
##############################
options(shiny.maxRequestSize=10000*1024^2)
source(paste(sourceDir, "/server.R", sep=""))
}

# Run the application
shinyApp(ui = ui, server = server)


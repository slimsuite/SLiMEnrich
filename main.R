#*********************************************************************************************************
#*********************************************************************************************************
# To run the app, please open slimenrichApp.R in RStudio and click "Run App".
# Further details can be found at the GitHub site.
################# ::: APP INFO ::: ######################
info = list(
  apptitle = "SLiMEnrich",
  version = "1.3.2",
  lastedit = "14 Jul 2018",
  author = "Sobia Idrees & Richard J. Edwards",
  contact = "richard.edwards@unsw.edu.au",
  description = "SLiMEnrich predicts Domain Motif Interactions (DMIs) from Protein-Protein Interaction (PPI) data and analyses enrichment through permutation test."
)
#*********************************************************************************************************
#*********************************************************************************************************
devmode = FALSE   # This affects some of the printing to screen
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
#V1.3.0 - Updated the DMI Strategy Handling. Added separate load functions. Modified progress bars and notifications.
#V1.3.1 - Partitioned out some general code into main.R. Updated ui.R and server.R. Updated slimenrich.R.
#       - Added options to slimenrich.R for output directory, DMI strategy, randomisations and histogram settings.
#V1.3.2 - Modified ui.R to load main.R for server functions. Tidied up library loading.
##############################
#SLiMEnrich program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# SLiMEnrich program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
# Package list. This has been consolidated for both the shiny and standalone implementations
package_names = c("shiny", "ggplot2", "colourpicker", "shinyBS", "shinythemes", "DT", "shinyjs", "visNetwork", "igraph","markdown","plotly", "plyr", "shinyWidgets")
if(devmode){
  load_or_install(package_names)
}else{
  suppressMessages(
    suppressWarnings(
      load_or_install(package_names)
  ))
}
##############################
#SETUP DATA
##############################
#i# This function is called by the sever (or server.R) in a reactiveValues() call.
#i# This way, it will be updated whenever the settings that affect it are changed.
#i# Data itself is stored in a list returned by the adata function, e.g. adata$data
#i# Individual elements are then accessed as adata$data$ELEMENT
setupData = function(){
  emptydb = list()
  for(rkey in c("PPI","Motifs","Domains","DMI")){
    emptydb[[rkey]] = data.frame()
  }
  emptydb$loads = list()   # Counter to register the Load button being pressed.
  for(rkey in c("PPI","Motifs","Domains","DMI","Calculate")){
    emptydb$loads[[rkey]] = 0
  }
  emptydb$activitylog = c(paste0("[",as.POSIXlt(Sys.time()),"] - ",info$apptitle," Version ",info$version,": Run started."))
  # Return data list
  return(emptydb)
}


####################################################################
### Functions for loading data
####################################################################
# Load PPI data (mProtein, dProtein)
#PPI2 <- loadPPIData(input)
loadPPIData <- function(input){
  PPIFile<-input$PPI
  if(is.null(PPIFile)){
    #showNotification("Either SLiM prediction or PPI file is missing. Loading example data", type = "warning", duration = NULL)
    PPI2<-read.csv("data/PPIs.csv",header=TRUE,sep=",")
  }
  else{
    
    PPI2<-read.csv(PPIFile$datapath,header=TRUE,sep=",")
  }
  return(PPI2)
}
parsePPIData <- function(input,PPI2){
  PPI2 <- PPI2[,c(input$ppimprotein,input$ppidprotein)]
  # Make unique pairs
  PPI2 = unique(PPI2)
  colnames(PPI2) = c("mProtein","dProtein")
  return(PPI2)
}

### Making the "Motif" table (mProtein-Motif)
#Motif <- loadDataMotif(input)
loadDataMotif <- function(input){
  #File upload check
  MotifFile<-input$Motif
  if(input$SLiMrunid){
    Motif<-read.delim(paste0("http://rest.slimsuite.unsw.edu.au/retrieve&jobid=",input$SLiMRun,"&rest=occ"),header=TRUE,sep=",")
  }
  else{
    # Check whether file loaded
    if(is.null(MotifFile)){
      fname <- "data/known.occ.csv"
    }else{
      fname <- MotifFile$datapath
    }
    # Check whether csv or tdt
    if(substr(fname,nchar(fname)-2,nchar(fname)) %in% c("csv","CSV")){
      Motif<-read.csv(fname,header=TRUE,sep=",")
    }else{
      Motif<-read.csv(fname,header=TRUE,sep="\t")
    }
  }
  return(Motif)
}
parseDataMotif <- function(input,Motif){
  if(input$SLiMrunid){
    motfield = "Motif"
    protfield = "AccNum"
  }
  else{
    # Select mProtein and Motif IDs
    # This will have input$motifmprotein and input$motifmotif text boxes 
    if(input$motifmprotein %in% colnames(Motif)){
      protfield = input$motifmprotein
    }else{
      protfield = "mProtein"
    }
    if(input$motifmotif %in% colnames(Motif)){
      motfield = input$motifmotif
    }else{
      motfield = "Motif"
    }
  }
  # NOTE: For direction dProtein links, this table will be replaced by duplicated DMI table fields
  # Pull out required columns dependent on strategy
  if(input$DMIStrategy %in% c("elmiprot")){
    # Direct protein links will use protein IDs from DMI file as domain IDs
    Motif <- loadDataMotifDomain(input)
    Motif <- parseDataMotifDomain(input,Motif)
    Motif <- Motif[,c("Motif","Motif")]
  }else{
    Motif <- Motif[,c(protfield,motfield)]
  }
  colnames(Motif) <- c("mProtein","Motif")
  return(Motif)
}


### Making the "Domain" table (Domain-dProtein)
#dProtein <- loadDatadomain(input)
loadDatadomain <- function(input){
  DomainFile<-input$domain
  # Check whether file loaded
  if(is.null(DomainFile)){
    fname <- "data/domain.csv"
  }else{
    fname <- DomainFile$datapath
  }
  # Check whether csv or tdt
  if(substr(fname,nchar(fname)-2,nchar(fname)) %in% c("csv","CSV")){
    dProtein<-read.csv(fname,header=TRUE,sep=",")
  }else{
    dProtein<-read.csv(fname,header=TRUE,sep="\t")
  }  
  return(dProtein)
}
parseDatadomain <- function(input,dProtein){
  # Select Domain and dProtein IDs
  # This will have input$domaindomain and input$domaindprotein text boxes 
  if(input$domaindomain %in% colnames(dProtein)){
    domfield = input$domaindomain
  }else{
    domfield = "Domain"
  }
  if(input$domaindprotein %in% colnames(dProtein)){
    protfield = input$domaindprotein
  }else{
    protfield = "dProtein"
  }
  
  # NOTE: For direction dProtein links, this table will be replaced by duplicated DMI table fields
  # Pull out required columns dependent on strategy
  if(input$DMIStrategy %in% c("elmiprot","elmcprot")){
    # Direct protein links will use protein IDs from DMI file as domain IDs
    dProtein <- loadDataMotifDomain(input)
    dProtein <- parseDataMotifDomain(input,dProtein)
    dProtein <- dProtein[,c("Domain","Domain")]
  }else{
    dProtein <- dProtein[,c(domfield,protfield)]
  }
  colnames(dProtein) <- c("Domain","dProtein")
  return(dProtein)
}


### Making the "DMI" table (Motif-Domain)
#Domain <- loadDataMotifDomain(input)
loadDataMotifDomain <- function(input){
  # Check whether data already loaded
  MotifDomainFile<-input$MotifDomain
  
  # Choose which file to use for DMI data
  # Default strategy is elmc-protein
  fname <- "data/elm_interactions.tsv"
  # Elm	Domain	interactorElm	interactorDomain	StartElm	StopElm	StartDomain	StopDomain	AffinityMin	AffinityMax	PMID	taxonomyElm	taxonomyDomain	
  # CLV_Separin_Fungi	PF03568	Q12158	Q03018	175	181	1171	1571	None	None	10403247,14585836	"559292"(Saccharomyces cerevisiae S288c)	"559292"(Saccharomyces cerevisiae S288c)
  if(input$DMIStrategy %in% c("elmcdom")){
    fname <- "data/motif-domain.tsv"
  }
  # "ELMidentifier"	"InteractionDomainId"	"Interaction Domain Description"	"Interaction Domain Name"
  #  "CLV_NRD_NRD_1"	"PF00675"	"Peptidase_M16"	"Insulinase (Peptidase family M16)"  
  # Check whether file loaded - over-rides default
  if(! is.null(MotifDomainFile)){
    fname <- MotifDomainFile$datapath
  }
  #writeLines(fname)
  
  # Check whether csv or tdt
  if(substr(fname,nchar(fname)-2,nchar(fname)) %in% c("csv","CSV")){
    Domain<-read.csv(fname,header=TRUE,sep=",")
  }else{
    Domain<-read.csv(fname,header=TRUE,sep="\t")
  }  
  #writeLines(paste(colnames(Domain)))
  return(Domain)
}
parseDataMotifDomain <- function(input,Domain){
  # Select Motif and Domain IDs
  # This will have input$dmimotif and input$dmidomain text boxes 
  if(input$dmimotif %in% colnames(Domain)){
    motfield = input$dmimotif
  }else{
    motfield = "Motif"
  }
  if(input$dmidomain %in% colnames(Domain)){
    domfield = input$dmidomain
  }else{
    domfield = "Domain"
  }
  Domain <- Domain[,c(motfield,domfield)]
  colnames(Domain) <- c("Motif","Domain")
  return(Domain)
}

#### Generate predicted DMis
# predDMIs <- generatePredictedDMIs(input)
generatePredictedDMIs <- function(input,appdata){
  
  Domain <- appdata$DMI #loadDataMotifDomain(input)
  dProtein <- appdata$Domains # loadDatadomain(input)
  Motif <- appdata$Motifs  #inputDataMotif()
  PPI2 <- appdata$PPI  # loadPPIData(input)
  
  Motif_NR<-unique(Motif)
  
  #Join/Merge two files based on Motif
  Join <- merge( Motif_NR, Domain, by="Motif")
  #print(Join)
  #joined both files based on Domains
  DMI <- merge(Join, dProtein,by="Domain")
  #Filtered unique DMIs
  Uni_DMI <- unique(DMI)
  #print(Uni_DMI)
  
  #PPI-DMI Mapping
  ########################################################################
  #names(PPI2) <- c("mProtein", "dProtein")
  predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
  Uni_predDMIs <- unique(predDMI)
  #names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
  predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
  #print(predDMIs)
  
  return(predDMIs)
}
##############################################################

##############################################################
# Simulating input settings for the commandline version
makeInputSettings <- function(){
  return(list(
    DMIStrategy = "elmcprot",
    PPI=list(datapath=paste0(rdir,"/data/PPIs.csv")),  # Default PPI file
    ppimprotein="mProtein",
    ppidprotein="dProtein",
    Motif=list(datapath=paste0(rdir,"/data/known.occ.csv")),
    SLiMrunid=FALSE,
    SLiMRun="",
    motifmprotein = "AccNum",
    motifmotif = "Motif",
    domain = list(datapath = paste0(rdir,"/data/domain.csv")),
    domaindomain = "pfam",
    domaindprotein = "accnum",
    MotifDomain = list(datapath = "ELM data"), # paste0(rdir,"/data/elm_interactions.tsv")),
    dmimotif = "Elm",
    dmidomain = "interactorDomain",
    shufflenum = 1000,
    xlimend = 0,
    binwidth = 1
  ))
}
##############################################################

##############################################################
# Documentation markdown
##############################################################
# This is a list containing vectors of markdown text. 
# These are joined with \n characters prior to rendering, e.g.
# each element in the vector is a line of markdown.
# These are used for pairs of htmlOutput() in ui.R and renderUI() in server.R.
# NOTE: An alternative is to load MD in the UI:  includeMarkdown("README.md")     
docs = list()
docs$tables = c("#### Uploaded Data","","Select PPI data, DMI strategy and optional data uploads.","Click **LOAD DATA**.",
                "Uploaded data will be displayed below.",
                "Check the **Show parsed data columns** box to display which data have been parsed from each input file.",
                "","**Note:** To analyse the example dataset, press **LOAD DATA** without uploading any files.",
                "","---","")

##############################################################

  







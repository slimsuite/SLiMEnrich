#!/usr/bin/env Rscript
#*********************************************************************************************************
#*********************************************************************************************************
# Please see main.R for App version, history and license information.
#*********************************************************************************************************
#*********************************************************************************************************
# Load main.R data and functions
#options(warn = -1)
# Check whether packages of interest are installed
is_installed = function(mypkg) is.element(mypkg, installed.packages()[,1])
if(!is_installed("stringr")) { install.packages("stringr") } 
library(stringr)
#devmode = FALSE   # This affects some of the printing to screen
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
rdir <- dirname(script.name)
rdir = str_replace_all(rdir,fixed("~+~")," ")
source(paste0(rdir,"/main.R"))
writeLines(c("","","########################################"))
writeLines(paste("### Running",info$apptitle,"Version",info$version,"###"))
writeLines(c("","","########################################")[3:1])
writeLines(c(paste("Run:", as.POSIXlt(Sys.time())),""))
#*********************************************************************************************************
#*********************************************************************************************************
# Additional packages for standalone version. Shiny server packages are loaded in main.R
package_names = c("optparse")
if(devmode){
  load_or_install(package_names)
}else{
  suppressMessages(
    suppressWarnings(
      load_or_install(package_names)
    ))
}
##########################################################################
## Options
input = makeInputSettings()
##########################################################################
# Commandline arguments
option_list = list(
  #PPI file
  make_option(c("-p", "--pFile"), type="character", default=NULL, 
              help="PPI file name [required]. Should have mProtein and dProtein fields.", metavar="FILENAME"),
  #SLiMs file 
  make_option(c("-s", "--mFile"), type="character", default=paste0(rdir,"/data/known.occ.csv"), 
              help="SLiM occurrence file name. Should have mProtein (or AccNum) and Motif fields. Not used if DMI strategy is `elmiprot`. [Default = ELM instances]", metavar="FILENAME"),
  #Motif-Domain file
  make_option(c("-m", "--mdFile"), type="character", default="ELM data", 
              help="DMI file name. Should have Motif and Domain fields unless using ELM input. [Default = ELM interaction data]", metavar="FILENAME"),
  #Domain file
  make_option(c("-d", "--dFile"), type="character", default=paste0(rdir,"/data/domain.csv"), 
              help="Domain file name. Should have Domain (or pfam) and dProtein (or accnum) fields. Not used if DMI strategy is `elmiprot` or `elmcprot`. [Default = Pfam domain for human Uniprot]", metavar="FILENAME"),
  #DMI strategy
  make_option(c("-t", "--strategy"), type="character", default="elmcprot", 
              help="DMI strategy (elmiprot/elmcprot/elmcdom). See docs for details.", metavar="character"),
  #Output directory
  make_option(c("-o", "--output"), type="character", default="./output/", 
              help="Output directory. Trailing characters will set file prefix (e.g. -o ./output/myrun.). Note: by default, the output directory is placed in the run directory. [Default = ./output/]", metavar="PATH"),
  #Randomisations
  make_option(c("-r", "--random"), type="integer", default=1000, 
              help="Number of PPI randomisations. [Default = 1000]", metavar="integer"),
  #Reuse data
  make_option(c("-u", "--userandom"), type="logical", default=FALSE, action = "store_true", 
              help="Whether to load and reuse existing PPI randomisations. WARNING: Make sure PPI input is unchanged when using this option. [Default = FALSE]", metavar=""),
  #X axis maximum
  make_option(c("-x", "--xmax"), type="integer", default=0, 
              help="Extend histogram x-axis to xmax.", metavar="integer"),
  #Histogram bin size
  make_option(c("-b", "--binsize"), type="integer", default=1, 
              help="Set histogram bin size. [Default = 1]", metavar="integer")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
#########################################################################
# Update input from options
if(opt$strategy %in% c("elmiprot","elmcprot","elmcdom")){
  input$DMIStrategy = opt$strategy
}else{
  writeLines(paste0("WARNING: DMI strategy (--strategy ",opt$strategy,") not recognised!"))
  writeLines(paste0("Using default: --strategy elmcprot"))
  input$DMIStrategy = "elmcprot"
}
# Set file names
#print(opt$mdFile)
if(opt$mdFile == "ELM data"){
  if(input$DMIStrategy %in% c("elmcdom")){
    input$MotifDomain$datapath = paste0(rdir,"/data/motif-domain.tsv")
    input$dmimotif = "ELMidentifier"
    input$dmidomain = "InteractionDomainId"  
  }else{
    input$MotifDomain$datapath = paste0(rdir,"/data/elm_interactions.tsv")
    if(input$DMIStrategy %in% c("elmiprot")){
      input$dmimotif = "interactorElm"
    }      
  }
}else{
  input$MotifDomain$datapath = opt$mdFile
}
input$PPI$datapath = opt$pFile
input$Motif$datapath = opt$mFile
input$domain$datapath = opt$dFile

#!# Add slimrun option SLiMrunid=FALSE,
#!# SLiMRun="",
# motifmprotein = "AccNum",
# motifmotif = "Motif",

input$shufflenum = opt$random

#!# Add checks for file existence

#*************************************************************************************
#Step 1: Potential DMIs
####################################################
adata = setupData()
writeLines('Loading data...')
# Load PPI data (mProtein-dProtein) -> adata$data$PPI
writeLines(paste0('Loading PPI from ', input$PPI$datapath))
adata$data$FullPPI = loadPPIData(input)
adata$data$PPI = unique(parsePPIData(input,adata$data$FullPPI))
if(devmode){ head(adata$data$PPI) }
# Making the "Motif" table (mProtein-Motif) -> adata$data$Motifs
writeLines(paste0('Loading Motifs from ', input$Motif$datapath))
adata$data$FullMotifs = loadDataMotif(input)
adata$data$Motifs = unique(parseDataMotif(input,adata$data$FullMotifs))
# Making the "DMI" table (Motif-Domain) -> adata$data$DMI
writeLines(paste0('Loading DMI from ', input$MotifDomain$datapath))
#print(input)
adata$data$FullDMI = loadDataMotifDomain(input)
adata$data$DMI = unique(parseDataMotifDomain(input,adata$data$FullDMI))
# Making the "Domain" table (Domain-dProtein) -> adata$data$Domains
writeLines(paste0('Loading domains from ', input$domain$datapath))
adata$data$FullDomains = loadDatadomain(input)
adata$data$Domains = unique(parseDatadomain(input,adata$data$FullDomains))
# Report loaded data
Data = adata$data
writeLines(paste0('PPI: ', nrow(Data$FullPPI), " (", ncol(Data$FullPPI), " fields); ", nrow(Data$PPI), " NR"))
writeLines(paste0('Motifs: ', nrow(Data$FullMotifs), " (", ncol(Data$FullMotifs), " fields); ", nrow(Data$Motifs), " NR"))
writeLines(paste0('DMI: ', nrow(Data$FullDMI), " (", ncol(Data$FullDMI), " fields); ", nrow(Data$DMI), " NR"))
writeLines(paste0('Domains: ', nrow(Data$FullDomains), " (", ncol(Data$FullDomains), " fields); ", nrow(Data$Domains), " NR"))

### Generate potential DMI
PPI2 <- Data$PPI
ppidProtein = as.character(unique(PPI2$dProtein))
ppimProtein = as.character(unique(PPI2$mProtein))
Domain <- Data$DMI
dProtein <- Data$Domains
dProtein <- dProtein[dProtein$dProtein %in% ppidProtein,]
Motif_NR <- Data$Motifs
Motif_NR <- Motif_NR[Motif_NR$mProtein %in% ppimProtein,]
#Join/Merge two files based on Motif
Join <- merge(Motif_NR, Domain, by="Motif")
#Domain-dProtein Mapping
DMI <- merge(Join, dProtein,by="Domain")
#Filtered unique DMIs
Uni_DMI <- unique(DMI)
#Named the header of output file
Uni_DMI <- Uni_DMI[, c("mProtein","Motif", "Domain", "dProtein")]
adata$data$potentialDMI = Uni_DMI
adata$data$potentialDMINR = unique(Uni_DMI[,c("mProtein","dProtein")])
writeLines(paste(nrow(adata$data$potentialDMI),"potential DMI;",nrow(adata$data$potentialDMINR),"NR"))

# Setup output directory
outdir = dirname(paste0(opt$output,"/"))
if(! file.exists(outdir) & outdir != "."){
  dir.create(outdir)
  writeLines(paste("A directory named",outdir, "has been created in the current working directory"))
}
output_potentialDMIs <- write.csv(Uni_DMI, paste0(opt$output,"potentialDMIs.csv"), row.names = FALSE)
writeLines("potentialDMIs.csv file has been saved in output folder")

#*************************************************************************************

#Step 2: Predicted DMIs
####################################################
#PPI-DMI Mapping                                                         
########################################################################

# PPI2<-read.csv(opt$pFile,header=TRUE,sep=",")
# PPI2 <- unique(PPI2)
# names(PPI2) <- c("mProtein", "dProtein")
writeLines("Finding predicted DMIs...", sep="\n")
predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
Uni_predDMIs <- unique(predDMI)
names(Uni_predDMIs) <- c("mProtein", "dProtein", "Motif", "Domain")
predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
if(devmode){ print(predDMIs) }
adata$data$predDMI = predDMIs
adata$data$predDMINR = unique(predDMIs[,c("mProtein","dProtein")])
writeLines(paste(nrow(adata$data$predDMI),"predicted DMI;",nrow(adata$data$predDMINR),"NR"))
output_predictedDMIs <- write.csv(predDMIs, paste0(opt$output,"predictedDMIs.csv"), row.names = FALSE)
writeLines("predictedDMIs.csv file has been saved in output folder")

#*************************************************************************************

#Step 3: Statistics
####################################################
writeLines("Calculating number of mProteins, motifs, domains and dProteins...", sep="\n")
png(filename=paste0(opt$output,"summaryStats.png"))
colors=c("cadetblue1", "deepskyblue2", "blue", "darkblue")
#Select unique Motifs
uniq_motif <- unique(predDMI$Motif)
a <- length(uniq_motif)
#Select unique Domain
uniq_domain <- unique(predDMI$Domain)
b <- length(uniq_domain)
#Select unique mProtein
uniq_mprotein <- unique(predDMI$mProtein)
c <- length(uniq_mprotein)
#Select unique dproteins
uniq_dprotein <- unique(predDMI$dProtein)
d <- length(uniq_dprotein)
uniq_count<-c(a = length(uniq_motif), b = length(uniq_domain), c = length(uniq_mprotein), d = length(uniq_dprotein))
#Create pie chart
x <- c(a,b,c,d)
#Label names for the chart
labels <- c("Motif", "Domain", "mProtein", "dProtein")
#Created bar char of the unique values
barplot(x, main="Statistics of DMIs", col = colors)
legend("topright", 
       legend = c(paste("Motifs=",a),paste("Domains=",b),paste("mProteins=",c),paste("dProteins=",d)), fill = colors)
noprint <- dev.off()
writeLines("Summary Statistics Bar chart has been saved in output folder")

     
#####################################################Enrichment Analysis####################################################
#*********************************************************************************************************************************
    
#*************************************************************************************
    
# Step 4: rPPI-DMI Mapping
####################################################
# This step generates 1000 random data and then compares each file with the potential DMIs dataset to generate a list of numbers (matches found in each PPI file).
# Randomized data will be stored in a new directory created in App folder named as "Randomdata" and can be accessed later if required.
# A file named randomNumbers will be generated in "Randomdata" that will be used to generate Histogram later.
####################################################

# Randomization/Permutations                                                  
      
#PPI_data<-read.csv(opt$pFile,header=TRUE,sep=",")
#PPI_data<-unique(PPI_data)
#names(PPI_data) <- c("mProtein", "dProtein")
PPI_data = adata$data$PPI
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

#rdir = paste0(opt$output,"Randomdata/")
rdir = paste0("Randomdata/",basename(opt$pFile),"/")

if(! file.exists(rdir)){
  dir.create(rdir)
  writeLines(paste("Created directory:",rdir))
}
writeLines(paste("Randomised PPI will be output to", rdir))
#for loop to create N randomized data
shufflenum = input$shufflenum
writeLines(paste("Generating",shufflenum,"randomised PPI datasets..."))
#!# Add re-use of randomises PPI based on file name
reuse = opt$userandom
reused = 0
for (j in 1:shufflenum) {
  jfile = paste0(rdir,"rPPI.",j,".csv")
  if(reuse & file.exists(jfile)){
    if(devmode){ writeLines(paste("Reusing",jfile))}
    reused = reused + 1
  }else{
    
    permutation<-PermutationFunction(PPI_Matrix, k = length(PPI_Matrix))
    permutation2<-PermutationFunction(PPI_Matrix2, k = length(PPI_Matrix2))
    final_file<- c(paste(permutation,permutation2, sep = ":"))
    newCol1<-strsplit(as.character(final_file),':',fixed=TRUE)
    df<-data.frame(final_file,do.call(rbind, newCol1))
    subset = df[,c(2,3)]
    names(subset)<-c("mProtein", "dProtein")
    write.csv(subset, jfile, row.names = FALSE)
  }
  if(j %% as.integer(shufflenum/10) == 0 & j < shufflenum){ writeLines(paste(j,"of",shufflenum,"...")) }
}
writeLines(paste0(shufflenum," shuffled PPI datasets have been created in Randomdata folder"))
if(reuse & reused > 0){
  writeLines(paste('NOTE: --userandom flag set;',reused,"of",shufflenum,"files were found and will be reused."))
}
if(reuse & reused < shufflenum){   # Some files made so need to repeat calculation
  reuse = FALSE
}
writeLines("Now mapping random PPI data to potential DMIs...")
#rPPI-DMI Mapping                                                               
#################################################################################
#!# Make devmode a toggle
#i# If devmode is on, will reuse calculated data
#!# Made this a devmode option because it is dangerous to reuse generic filenames when data can be different
randfile = paste0(opt$output,"randomDMI.",input$DMIStrategy,".csv")
if(reuse & file.exists(randfile)){
  x<-read.csv(randfile, sep = ",", header = FALSE)
  writeLines("NOTE: --userandom flag set and randomNumbers file already exists. Will load and reuse. Ensure input is unchanged.")
}else {
  for (i in 1:shufflenum) {
    rPPI <- read.table(paste0(rdir,"rPPI.", i, ".csv"),
                       stringsAsFactors=FALSE, sep=",", strip.white=TRUE)
    names(rPPI)<-c("mProtein", "dProtein")
    names(Uni_DMI) <- c("mProtein","Motif", "Domain","dProtein")
    DMI_rPPI <- merge(Uni_DMI, rPPI, by= c("mProtein", "dProtein"))
    DMI_rPPI <- unique(DMI_rPPI[,c("mProtein", "dProtein")])
    Matches <- nrow(DMI_rPPI)
    if(devmode){
      writeLines(paste0("File ",i,": ",Matches))
    }else{
      if(i %% as.integer(shufflenum/10) == 0 & i < shufflenum){ writeLines(paste(i,"of",shufflenum,"...")) }
    }
    output_randomNumbers <- write.table(Matches, randfile, col.names = FALSE, append = i > 1, row.names = FALSE)
  }
}

#*************************************************************************************
    
#Step 5: Histogram
####################################################
# Plotting/comparing NR predicted DMIs
predDMIs = adata$data$predDMINR
predDMInum = nrow(predDMIs)
input$binwidth = opt$binsize
input$xlimend = opt$xmax

png(filename=paste0(opt$output,"Histogram.png"), width = 1200, height = 800, units = "px", pointsize = 12)
x<-read.csv(randfile, sep = ",", header = FALSE)
names(x) <- "values"
#bins<- seq(min(x), max(x), length.out = 20 + 1)
xlimmax = max(c(input$xlimend,max(x)*1.1,predDMInum*1.1))
bins<- seq(0, xlimmax+input$binwidth, input$binwidth)
ymax = max(hist(x$values, breaks=bins, plot=FALSE)$counts) * 1.5   #shufflenum/2

h<- hist(x$values, breaks=bins, col = "red", border = 'black', main="Distribution of Random DMIs", ylab="Frequency of Random DMIs",
         xlab="Number of Random DMIs", ylim = c(0,ymax), xlim = c(0, xlimmax), cex.main=1.5, cex.lab=1.5,cex.axis=1.5)

xfit <- seq(min(x$values),max(x$values),length=40)
yfit <- dnorm(xfit, mean=mean(x$values),sd=sd(x$values))
yfit<- yfit*diff(h$mids[1:2])*length(x$values)
lines(xfit,yfit)

#axis(side=3, lwd = 0, lwd.ticks = 4, at=nrow(predictedDMIs()), lend=1, labels = FALSE, tcl=5, font=2, col = "black", padj = 0, lty = 3)
#shows the observed value
mtext(paste("Observed value is: ", predDMInum), side = 3, at=nrow(predDMIs), font = 4)

pvalue = length(x[x >= predDMInum])/shufflenum
pplace = xlimmax
if(predDMInum > xlimmax/2){ pplace = predDMInum/2 }

if(pvalue == 0){ 
  mtext(paste0("P-value is: < ", 1/shufflenum), side = 3, at=pplace, font = 4, col = "red")
  pvalue <- paste0("< ", 1/input$shufflenum)
}else{
  #mtext(paste0("P-value is: ", pvalue), side = 3, at=predDMInum+90, font = 4, col = "red")
  mtext(paste0("P-value is: ", pvalue), side = 3, at=pplace, font = 4, col = "red")
}

#points arrow on the observed value
arrows(predDMInum, ymax, predDMInum, 0, lwd = 2, col = "black", length = 0.1, lty = 3)
noprint <- dev.off()
ObsDMI <- nrow(adata$data$predDMI)
ObsDMINR <- predDMInum
MeanRandDMINR <- mean(x$values)
#pvalue <-   length(x[x >= predDMInum])/shufflenum
FDR <- mean(x$values)/nrow(predDMIs)
Escore <- nrow(predDMIs)/mean(x$values)
#!# Improve this: transpose and add file info
sum <- data.frame(ObsDMI, ObsDMINR, MeanRandDMINR, pvalue, FDR, Escore)
summary <- write.csv(sum, paste0(opt$output,"summary.csv"), row.names = FALSE, quote = FALSE)
head(sum)
writeLines(c("","Summary table saved to output directory"))

#!#writeLines("summary Table has been created in output directory")
#*************************************************************************************
    
#Step 6: DMI Network
####################################################
predDMIs = adata$data$predDMI
#Network
first <- predDMIs[,c("mProtein","Motif")]
g <- graph.data.frame(first, directed = F)
#igraph.options(plot.layout=layout.graphopt, vertex.size=10)
V(g)$color <- ifelse(V(g)$name %in% predDMI[,1], "#A93226", "#F7DC6F")
V(g)$shape <- ifelse(V(g)$name %in% predDMI[,1], "square", "box")
#plot(g,  edge.color="orange")
#visIgraph(g)
##print(first)
#V(g)$color <- "red"
#Second
second <- predDMIs[,c("Motif","Domain")]
g2 <- graph.data.frame(second, directed = F)
#igraph.options(plot.layout=layout.graphopt, vertex.size=10)
V(g2)$color <- ifelse(V(g2)$name %in% predDMI[,1], "#9B59B6", "#D35400")
V(g2)$shape <- ifelse(V(g2)$name %in% predDMI[,1], "box", "circle")
#plot(g2,  edge.color="orange")
#print(second)
#visIgraph(g2)
#V(g2)$color <- "green"
#Third

third <- predDMIs[,c("Domain","dProtein")]
g3 <- graph.data.frame(third, directed = F)
#igraph.options(plot.layout=layout.graphopt, vertex.size=10)
V(g3)$color <- ifelse(V(g3)$name %in% predDMI[,3], "#9B59B6", "#85C1E9")
V(g3)$shape <- ifelse(V(g3)$name %in% predDMI[,3], "vrectangle", "circle")

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
g6 <- visIgraph(g5,layout = "layout_nicely", physics = FALSE, smooth = TRUE, type = "square")

#visSave(g6, file = paste0(opt$output,"network.html"), selfcontained = TRUE)
visSave(g6, file = "network.html", selfcontained = TRUE)
noprint <- file.rename("network.html",paste0(opt$output,"network.html"))
writeLines("A DMI network has been saved as html file")
        
writeLines(c("",paste0("End ",info$apptitle," V",info$version," Run: ", as.POSIXlt(Sys.time()))))
        

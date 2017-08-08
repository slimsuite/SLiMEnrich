#*********************************************************************************************************
#*********************************************************************************************************
# Short Linear Motif Enrichment Analysis App (SLiMEnrich)
# Developer: **Sobia Idrees**
# Version: 1.0.2
# Description: SLiMEnrich predicts Domain Motif Interactions (DMIs) from Protein-Protein Interaction (PPI) data and analyzes enrichment through permutation test.
#*********************************************************************************************************
#*********************************************************************************************************
##############################
#Version History
##############################
#V1.0.1 - Generic naming
#V1.0.2-  Added commandline arguments to select files
##############################
# Argument Description
#SLiMs file option
# -s
#Domain file option
# -d
#Motif-Domain file option
# -m
#PPI file argument
# -p
##############################
#Required Libraries
##############################
#!/usr/bin/env Rscript
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
load_or_install(c("ggplot2", "visNetwork", "igraph","optparse"))
##########################################################################
## Options
# Commandline arguments
##########################################################################
option_list = list(
  #SLiMs file 
  make_option(c("-s", "--mfile"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  #Domain file
  make_option(c("-d", "--dFile"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  #Motif-Domain file
  make_option(c("-m", "--mdFile"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  #PPI file
  make_option(c("-p", "--pFile"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#########################################################################

#*************************************************************************************
#Step 1: Potential DMIs
####################################################
#Read uploaded data
Motif<-read.csv(opt$mfile,header=TRUE,sep=",")[,c('AccNum','Motif')]
Motif_NR<-unique(Motif)
#If motif-domain file is not uploaded then read from data folder
if(is.null(opt$mdFile)){
  Domain<-read.csv("data/motif-domain.tsv",header=TRUE,sep="\t")[,c(1:2)]
}else{
  Domain<-read.csv(opt$mdFile,header=TRUE,sep="\t")[,c(1:2)]
}
#If domain file is not uploaded then read from data folder
if(is.null(opt$dFile)){
  dProtein<-read.csv("data/domain.csv",header=TRUE,sep=",")[,c('pfam','accnum')]
}else{
  dProtein<-read.csv(opt$dFile,header=TRUE,sep=",")[,c('pfam','accnum')]
}
#Motif-Domain Mapping                                            
#Rename the columns in two data
names(Motif_NR) <- c("Seq", "Motif")
names(Domain) <- c("Motif", "Domain")
#Join/Merge two data based on Motif                 
Join <- merge(Motif_NR, Domain, by="Motif")
#print(Join)
names(Join) <- c("Motif", "Seq", "Domain")  #Change header of the output file
#Domain-dProtein Mapping                                       
#Load results from the previous code)
names(dProtein) <- c("Domain", "dProteins")
writeLines("Finding potential DMIs...", sep="\n")
#joined both data based on Domain
DMI <- merge(Join, dProtein,by="Domain")
#Filtered unique DMIs
Uni_DMI <- unique(DMI)
#Named the header of output file
names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
Uni_DMI <- Uni_DMI[, c("mProtein","Motif", "Domain", "dProtein")]
print(Uni_DMI)
dir.create("output")
print("A directory named output has been created in the current working directory")
output_potentialDMIs <- write.csv(Uni_DMI, "output/potentialDMIs.csv", row.names = FALSE)
print("potentialDMIs.csv file has been saved in output folder")

#*************************************************************************************

#Step 2: Predicted DMIs
####################################################
#PPI-DMI Mapping                                                         
########################################################################

PPI2<-read.csv(opt$pFile,header=TRUE,sep=",")
names(PPI2) <- c("mProtein", "dProtein")
writeLines("Finding predicted DMIs...", sep="\n")
predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
Uni_predDMIs <- unique(predDMI)
names(Uni_predDMIs) <- c("mProtein", "dProtein", "Motif", "Domain")
predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
print(predDMIs)    
output_predictedDMIs <- write.csv(predDMIs, "output/predictedDMIs.csv", row.names = FALSE)
print("predictedDMIs.csv file has been saved in output folder")

#*************************************************************************************

#Step 3: Statistics
####################################################
writeLines("Calculating number of mProteins, motifs, domains and dProteins...", sep="\n")
png(filename="output/summaryStats.png")
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
print("Summary Statistics Bar chart has been saved in output folder")
dev.off()

     
    #####################################################Enrichment Analysis####################################################
    #*********************************************************************************************************************************
    
    
    #*************************************************************************************
    
    #Step 4: rPPI-DMI Mapping
    ####################################################
    #This step generates 1000 random data and then compares each file with the potential DMIs dataset to generate a list of numbers (matches found in each PPI file).
    #Randomized data will be stored in a new directory created in App folder named as "Randomdata" and can be accessed later if required.
    #A file named randomNumbers will be generated in "Randomdata" that will be used to generate Histogram later.
    ####################################################

      #Randomization/Permutations                                                  
      
      PPI_data<-read.csv(opt$pFile,header=TRUE,sep=",")
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
      
      
      dir.create("Randomdata")
      
     # dirName <- paste0("Randomdata_", strsplit(as.character(PPIFile$name), '.csv'))
    #  dir.create(dirName)
      #for loop to create 1000 randomized data
      for (j in 1:1000) {
        permutation<-PermutationFunction(PPI_Matrix, k = length(PPI_Matrix))
        permutation2<-PermutationFunction(PPI_Matrix2, k = length(PPI_Matrix2))
        final_file<- c(paste(permutation,permutation2, sep = ":"))
        newCol1<-strsplit(as.character(final_file),':',fixed=TRUE)
        df<-data.frame(final_file,do.call(rbind, newCol1))
        subset = df[,c(2,3)]
        names(subset)<-c("mProtein", "dProtein")
        write.csv(subset, paste0("Randomdata/rPPI",j,".csv"), row.names = FALSE)
      }
      print("1000 Randomdata have been created in Randomdata folder")
      print("Now mapping random PPI data to potential DMIs")
      #rPPI-DMI Mapping                                                               
      #################################################################################
      if(file.exists("Randomdata/randomNumbers.csv")){
        x<-read.csv("Randomdata/randomNumbers.csv", sep = ",", header = FALSE)
        print("randomNumbers file already exists. Will load instead")
      }else {
      for (i in 1:1000) {
        rPPI <- read.table(paste0("Randomdata/rPPI", i, ".csv"),
                           stringsAsFactors=FALSE, sep=",", strip.white=TRUE)
        names(rPPI)<-c("mProtein", "dProtein")
        names(Uni_DMI) <- c("mProtein","Motif", "Domain","dProtein")
        DMI_rPPI <- merge(Uni_DMI, rPPI, by= c("mProtein", "dProtein"))
        Matches <- nrow(DMI_rPPI)
        print(paste0("File ",i,": ",Matches))
        
        output_randomNumbers <- write.table(Matches, "Randomdata/randomNumbers.csv", col.names = FALSE, append = TRUE, row.names = FALSE)
        
      }
        }
    #*************************************************************************************
    
    #Step 5: Histogram
    ####################################################
      png(filename="output/Histogram.png", width = 1200, height = 800, units = "px", pointsize = 12)
      x<-read.csv("Randomdata/randomNumbers.csv", sep = ",", header = FALSE)
      names(x) <- "values"
      bins<- seq(min(x), max(x), length.out = 20 + 1)
      h<- hist(x$values, breaks=bins, col = "red", border = 'black', main="Distribution of Random DMIs", ylab="Frequency of Random DMIs",
               xlab="Number of Random DMIs", ylim = c(0,500), xlim = c(0, max(x)+100), cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
      xfit <- seq(min(x$values),max(x$values),length=40)
      
      yfit <- dnorm(xfit, mean=mean(x$values),sd=sd(x$values))
      
      
      yfit<- yfit*diff(h$mids[1:2])*length(x$values)
      
      
      lines(xfit,yfit)
      #axis(side=3, lwd = 0, lwd.ticks = 4, at=nrow(predictedDMIs()), lend=1, labels = FALSE, tcl=5, font=2, col = "black", padj = 0, lty = 3)
      #shows the observed value
      mtext(paste("Observed value is: ", nrow(predDMIs)), side = 3, at=nrow(predDMIs), font = 4)
      mtext(paste0("P-value is: ", length(x[x >= nrow(predDMIs)])/1000), side = 3, at=nrow(predDMIs)+90, font = 4, col = "red")
      #points arrow on the observed value
      arrows(nrow(predDMIs), 500, nrow(predDMIs), 0, lwd = 2, col = "black", length = 0.1, lty = 3)
      dev.off()
      pvalue <-   length(x[x >= nrow(predDMIs)])/1000
      meanvalue <- mean(x$values)
      FDR <- mean(x$values)/nrow(predDMIs)
     sum <- data.frame(pvalue, meanvalue, FDR)
     summary <- write.csv(sum, "output/summary.csv", row.names = FALSE)
     print("summary Table has been created in output directory")
    #*************************************************************************************
    
    #Step 6: DMI Network
    ####################################################
   
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
        
        visSave(g6, file = "network.html", selfcontained = FALSE)
        print("A DMI network has been saved as html file")

#*********************************************************************************************************
#*********************************************************************************************************
# Short Linear Motif Enrichment Analysis App (SLiMEnrich)
# Developer: **Sobia Idrees**
# Version: 1.0.0
# Description: SLiMEnrich predicts Domain Motif Interactions (DMIs) from Protein-Protein Interaction (PPI) data and analyzes enrichment through permutation test.
#*********************************************************************************************************
#*********************************************************************************************************
##############################
#Required Libraries
##############################
library(ggplot2)
library(visNetwork)
library(igraph)
#########################################################################

#*************************************************************************************
#Step 1: Potential DMIs
####################################################
      #Read uploaded files
      ELM<-read.csv("Files/slimprob.occ.csv",header=TRUE,sep=",")[,c('AccNum','Motif')]
      ELM_NR<-unique(ELM)
      Pfam<-read.csv("Files/ELMs_Pfams.tsv",header=TRUE,sep="\t")[,c(1:2)]
      hProtein<-read.csv("Files/Pfams.csv",header=TRUE,sep=",")[,c('pfam','accnum')]
      #hProtein<-read.csv("Files/hProteins.csv",header=TRUE,sep=",")
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
      #print(Uni_DMI)
      dir.create("output")
      output_potentialDMIs <- write.csv(Uni_DMI, "output/potentialDMIs.csv", row.names = FALSE)
      print("potentialDMIs.csv file has been saved in output folder")
    #*************************************************************************************
    
    #Step 2: Predicted DMIs
    ####################################################
      #vhPPI-DMI Mapping                                                         
      ########################################################################
      PPI2<-read.csv("Files/PPIs.csv",header=TRUE,sep=",")
      names(PPI2) <- c("Protein1", "Protein2")
      predDMI <- merge(PPI2, Uni_DMI, by= c("Protein1", "Protein2"))
      Uni_predDMIs <- unique(predDMI)
      names(Uni_predDMIs) <- c("Protein1", "Protein2", "ELM", "Pfam")
      predDMIs <- Uni_predDMIs[, c("Protein1","ELM", "Pfam", "Protein2")]
      #print(predDMIs)    
      output_predictedDMIs <- write.csv(predDMIs, "output/predictedDMIs.csv", row.names = FALSE)
      print("predictedDMIs.csv file has been saved in output folder")
    
    #*************************************************************************************
    
    #Step 3: Statistics
    ####################################################
        
        png(filename="output/summaryStats.png")
        colors=c("cadetblue1", "deepskyblue2", "blue", "darkblue")
        #Select unique ELM
        uniq_elm <- unique(predDMI$ELM)
        a <- length(uniq_elm)
        #Select unique Pfam
        uniq_pfam <- unique(predDMI$Pfam)
        b <- length(uniq_pfam)
        #Select unique vORF
        uniq_vorf <- unique(predDMI$Protein1)
        c <- length(uniq_vorf)
        #Select unique hproteins
        uniq_hprotein <- unique(predDMI$Protein2)
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
        print("Summary Statistics Bar chart has been saved in output folder")
        dev.off()
     
    #####################################################Enrichment Analysis####################################################
    #*********************************************************************************************************************************
    
    
    #*************************************************************************************
    
    #Step 4: rPPI-DMI Mapping
    ####################################################
    #This step generates 1000 random files and then compares each file with the potential DMIs dataset to generate a list of numbers (matches found in each PPI file).
    #Randomized files will be stored in a new directory created in App folder named as "RandomFiles" and can be accessed later if required.
    #A file named randomNumbers will be generated in "RandomFiles" that will be used to generate Histogram later.
    ####################################################

      #Randomization/Permutations                                                  
      
      PPI_data<-read.csv("Files/PPIs.csv",header=TRUE,sep=",")
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
      
      
      dir.create("RandomFiles")
      
     # dirName <- paste0("RandomFiles_", strsplit(as.character(PPIFile$name), '.csv'))
    #  dir.create(dirName)
      #for loop to create 1000 randomized files
      for (j in 1:1000) {
        permutation<-PermutationFunction(PPI_Matrix, k = length(PPI_Matrix))
        permutation2<-PermutationFunction(PPI_Matrix2, k = length(PPI_Matrix2))
        final_file<- c(paste(permutation,permutation2, sep = ":"))
        newCol1<-strsplit(as.character(final_file),':',fixed=TRUE)
        df<-data.frame(final_file,do.call(rbind, newCol1))
        subset = df[,c(2,3)]
        names(subset)<-c("vProtein", "hProtein")
        write.csv(subset, paste0("RandomFiles/rPPI",j,".csv"), row.names = FALSE)
      }
      print("1000 RandomFiles have been created in RandomFiles folder")
      print("Now mapping random PPI files to potential DMIs")
      #rPPI-DMI Mapping                                                               
      #################################################################################
      if(file.exists("RandomFiles/randomNumbers.csv")){
        x<-read.csv("RandomFiles/randomNumbers.csv", sep = ",", header = FALSE)
        print("randomNumbers file already exists. Will load instead")
      }else {
      for (i in 1:1000) {
        rPPI <- read.table(paste0("RandomFiles/rPPI", i, ".csv"),
                           stringsAsFactors=FALSE, sep=",", strip.white=TRUE)
        names(rPPI)<-c("vProtein", "hProtein")
        names(Uni_DMI) <- c("vProtein","ELM", "Pfam","hProtein")
        DMI_rPPI <- merge(Uni_DMI, rPPI, by= c("vProtein", "hProtein"))
        Matches <- nrow(DMI_rPPI)
        print(paste0("File ",i,": ",Matches))
        
        output_randomNumbers <- write.table(Matches, "RandomFiles/randomNumbers.csv", col.names = FALSE, append = TRUE, row.names = FALSE)
        
      }
        }
    #*************************************************************************************
    
    #Step 5: Histogram
    ####################################################
      png(filename="output/Histogram.png", width = 1200, height = 800, units = "px", pointsize = 12)
      x<-read.csv("RandomFiles/randomNumbers.csv", sep = ",", header = FALSE)
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
        first <- predDMIs[,c("Protein1","ELM")]
        g <- graph.data.frame(first, directed = F)
        #igraph.options(plot.layout=layout.graphopt, vertex.size=10)
        V(g)$color <- ifelse(V(g)$name %in% predDMI[,1], "#A93226", "#F7DC6F")
        V(g)$shape <- ifelse(V(g)$name %in% predDMI[,1], "square", "box")
        #plot(g,  edge.color="orange")
        #visIgraph(g)
        ##print(first)
        #V(g)$color <- "red"
        #Second
        second <- predDMIs[,c("ELM","Pfam")]
        g2 <- graph.data.frame(second, directed = F)
        #igraph.options(plot.layout=layout.graphopt, vertex.size=10)
        V(g2)$color <- ifelse(V(g2)$name %in% predDMI[,1], "#9B59B6", "#D35400")
        V(g2)$shape <- ifelse(V(g2)$name %in% predDMI[,1], "box", "circle")
        #plot(g2,  edge.color="orange")
        #print(second)
        #visIgraph(g2)
        #V(g2)$color <- "green"
        #Third
        
        third <- predDMIs[,c("Pfam","Protein2")]
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

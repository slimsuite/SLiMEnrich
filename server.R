#*********************************************************************************************************
#*********************************************************************************************************
# Please see main.R for App version, history and license information.
#devmode = FALSE   # This affects some of the printing to screen
#thisisshiny = TRUE  # This determines the packages loaded.
source("main.R")
package_names = c("shiny", "ggplot2", "colourpicker", "shinyBS", "shinythemes", "DT", "shinyjs", "visNetwork", "igraph","markdown","plotly", "plyr", "shinyWidgets")
# if(devmode){
#   load_or_install(package_names)
# }else{
#   suppressMessages(
#     suppressWarnings(
#       load_or_install(package_names)
#     ))
# }
writeLines(paste("Running",info$apptitle,"Version",info$version))
#*********************************************************************************************************
#*********************************************************************************************************
#*********************************************************************************************************
#*********************************************************************************************************

##############################
# Server logic
##############################
options(shiny.maxRequestSize=10000*1024^2)
server <- shinyServer(function(input, output, session){
  # Setup the data list that will store loaded data
  adata <- reactiveValues(
    data = setupData()
  )
  # This function computes a new data set. It can optionally take a function,
  # updateProgress, which will be called as each row of data is added.
  #?# What is this for in the app?!
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
    toggle(id = "advsettings", anim = TRUE)
  })
  observeEvent(input$SLiMrunid, {
    toggle(id = "slimrun", anim = TRUE)
    if(input$SLiMrunid){
      hide("slimf")
    }
    else{
      show("slimf")
    }
  })
  
  
  #lowercase function
  lowername <- function(x) {
    colnames(x) <- tolower(colnames(x))
    x
    
  }
  # Generate notifications of loaded data; calculates predicted and potential DMI for later display
  observeEvent(input$run, {
    withProgress(message = 'Loading data', detail = "PPI", value = 0, {
      PPIFile<-input$PPI
      incProgress(0.25, detail = "DMI")
      MotifDomainFile<-input$MotifDomain
      incProgress(0.25, detail = "Motifs")
      MotifFile<-input$Motif
      incProgress(0.25, detail = "Domains")
      DomainFile<-input$domain
      incProgress(0.25, detail = "Complete")
    })
    
    if(is.null(PPIFile)){
      showNotification("PPI file not provided: loading example Adenovirus PHISTO dataset", type = "warning", duration = 5)
    }
    SliMJobId <- input$SLiMRun
    if(is.null(MotifFile) && SliMJobId == "" ){
      showNotification("SLiM file not provided: loading ELM instances", type = "warning", duration = 5)
    }
    if(is.null(MotifDomainFile)){
      showNotification("DMI file not provided: loading ELM interaction data", type = "warning", duration = 5)
    }
    if(is.null(DomainFile)){
      showNotification("Domain file not provided: loading reviewed human Uniprot Pfam data", type = "warning", duration = 5)
    }
    if(is.null(PPIFile)){
      showNotification("Loaded example Adenovirus PHISTO dataset", type = "warning", duration = 30)    
    }
    # Generate potential and predicted DMIs on loading data
    withProgress(message = 'Calculating DMI', detail = "Potential", value = 0, {
      potentialDMIs()
      showNotification(paste(nrow(adata$data$potentialDMI),"potential DMI;",nrow(adata$data$potentialDMINR),"NR"), type = "message", duration = NULL)
      incProgress(0.5, detail = "Predicted")
      predictedDMIs()
      showNotification(paste(nrow(adata$data$predDMI),"predicted DMI;",nrow(adata$data$predDMINR),"NR"), type = "message", duration = NULL)
      incProgress(0.5, detail = "Complete")
    })
    adata$data$loads$Calculate = input$run
  })
  
  # observeEvent(input$DMIStrategy) to update default DMI fields
  observeEvent(input$DMIStrategy, {
    #textInput(inputId="dmimotif",label = "DMI file Motif column", value = "ELM"),
    #textInput(inputId="dmidomain",label = "DMI file Domain column", value = "interactorDomain"),
    if(input$DMIStrategy == "elmiprot"){   # Motif=mProtein, Domain=dProtein
      updateTextInput(session, "dmimotif",
                      label = "DMI file Motif (mProtein) column",
                      value = "interactorElm"
      )
      updateTextInput(session, "dmidomain",
                      label = "DMI file Domain (dProtein) column",
                      value = "interactorDomain"
      )
    }
    if(input$DMIStrategy == "elmcprot"){   # Motif=mProtein, Domain=dProtein
      updateTextInput(session, "dmimotif",
                      label = "DMI file Motif column",
                      value = "Elm"
      )
      updateTextInput(session, "dmidomain",
                      label = "DMI file Domain (dProtein) column",
                      value = "interactorDomain"
      )
    }
    if(input$DMIStrategy == "elmcdom"){   # Motif=mProtein, Domain=dProtein
      updateTextInput(session, "dmimotif",
                      label = "DMI file Motif column",
                      value = "ELMidentifier"
      )
      updateTextInput(session, "dmidomain",
                      label = "DMI file Domain column",
                      value = "InteractionDomainId"
      )
    }
  })
  
  #####################################################Domain-Motif Interactions####################################################
  #*********************************************************************************************************************************
  #Uploaded Data
  ####################################################
  
  #i# DMI Strategies: input$DMIStrategy
  # "Link mProteins directly to dProteins (ELMi-Protein)" = "elmiprot", 
  # "Link Motif classes directly to dProteins (ELMc-Protein)" = "elmcprot",
  # "Link Motif classes to binding domains (ELMc-Domain)" = "elmcdom"),
  
  inputDataPPI <-eventReactive(input$run, {
    #return(loadPPIData(input))
    if(adata$data$loads$PPI < input$run){
      adata$data$loads$PPI = input$run
      adata$data$FullPPI = loadPPIData(input)
      adata$data$PPI = unique(parsePPIData(input,adata$data$FullPPI))
    }
    return(adata$data$PPI)
  })
  
  ### Making the "Motif" table (mProtein-Motif)
  inputDataMotif <-eventReactive(input$run, {
    #return(loadDataMotif(input))
    if(adata$data$loads$Motifs < input$run){
      adata$data$loads$Motifs = input$run
      adata$data$FullMotifs = loadDataMotif(input)
      adata$data$Motifs = unique(parseDataMotif(input,adata$data$FullMotifs))
    }
    return(adata$data$Motifs)
  })
  
  ### Making the "Domain" table (Domain-dProtein)
  inputDatadomain <-eventReactive(input$run, {
    #return(loadDatadomain(input))
    if(adata$data$loads$Domains < input$run){
      adata$data$loads$Domains = input$run
      adata$data$FullDomains = loadDatadomain(input)
      adata$data$Domains = unique(parseDatadomain(input,adata$data$FullDomains))
    }
    return(adata$data$Domains)
  })
  
  ### Making the "DMI" table (Motif-Domain)
  inputDataMotifDomain <-eventReactive(input$run, {
    if(adata$data$loads$DMI < input$run){
      adata$data$loads$DMI = input$run
      adata$data$FullDMI = loadDataMotifDomain(input)
      adata$data$DMI = unique(parseDataMotifDomain(input,adata$data$FullDMI))
    }
    return(adata$data$DMI)
  })
  
  #########################################################
  #Displaying Data Tables
  ####################################################
  # PPI Table display
  output$udata2<-renderDataTable({
    inputDataPPI()
    if(input$parseddata){
      adata$data$PPI
    }else{
      adata$data$FullPPI
    }
  },
  caption = tags$h4(tags$strong("Interaction File"))
  )
  # DMI Table display
  output$udata4<-renderDataTable({
    inputDataMotifDomain()
    if(input$parseddata){
      adata$data$DMI
    }else{
      adata$data$FullDMI
    }
  },
  caption = tags$h4(tags$strong("Motif-Domain File"))
  )
  # Motifs Table display 
  output$udata<-renderDataTable({
    inputDataMotif()
    if(input$parseddata){
      adata$data$Motifs
    }else{
      adata$data$FullMotifs
    }
  },
  caption = tags$h4(tags$strong("SLiM File"))
  )
  # Domains Table display 
  output$udata3<-renderDataTable({
    inputDatadomain()
    if(input$parseddata){
      adata$data$Domains
    }else{
      adata$data$FullDomains
    }
  },
  caption = tags$h4(tags$strong("Domain File"))
  )
  #*************************************************************************************
  #Step 1: Potential DMIs
  ####################################################
  potentialDMIs <-eventReactive(input$run, {
    #i# DMI Strategies: input$DMIStrategy
    # "Link mProteins directly to dProteins (ELMi-Protein)" = "elmiprot", 
    # "Link Motif classes directly to dProteins (ELMc-Protein)" = "elmcprot",
    # "Link Motif classes to binding domains (ELMc-Domain)" = "elmcdom"),
    # NOTE: The inputData functions now return unique two-field data tables
    writeLines("Potential DMI")
    if(adata$data$loads$Calculate < input$run){
      withProgress(message = 'Potential DMI', detail = "(DMI)", value = 0, {
        Domain <- inputDataMotifDomain()
        incProgress(0.2, detail = "(Domains)")
        dProtein <- inputDatadomain()
        incProgress(0.2, detail = "(Motifs)")
        Motif_NR <- inputDataMotif()
        #Join/Merge two files based on Motif
        incProgress(0.2, detail = "(Joining data)")
        Join <- merge(Motif_NR, Domain, by="Motif")
        #Domain-dProtein Mapping
        #joined both files based on domain
        incProgress(0.2, detail = "(Joining data)")
        DMI <- merge(Join, dProtein,by="Domain")
        #Filtered unique DMIs
        incProgress(0.2, detail = "(Filtering data)")
        Uni_DMI <- unique(DMI)
        #Named the header of output file
        Uni_DMI <- Uni_DMI[, c("mProtein","Motif", "Domain", "dProtein")]
        print(head(Uni_DMI))
        adata$data$potentialDMI = Uni_DMI
        adata$data$potentialDMINR = unique(Uni_DMI[,c("mProtein","dProtein")])
      })
    }
    return(adata$data$potentialDMI)
  })
  
  
  #shows the data table
  output$data<-DT::renderDataTable({
    #Run it only when run button is active
    input$nrpotdmi
    if(input$run)
    {
      #Generates progress bar
      style <- isolate(input$style)
      
      if(input$nrpotdmi){
        potDMI <- unique(potentialDMIs()[,c("mProtein","dProtein")])
        print(head(potDMI))
        formatStyle(datatable(potDMI), columns = 1:2, color = "black")
      }else{
        formatStyle(datatable(potentialDMIs()), columns = 1:4, color = "black")
      }
    }
    
  })
  #*************************************************************************************
  
  #Step 2: Predicted DMIs
  ####################################################
  predictedDMIs <-eventReactive(input$run, {
    writeLines("Predicting DMIs")
    if(adata$data$loads$Calculate < input$run){
      withProgress(message = 'Predicted DMI', detail = "(PPI)", value = 0, {
        PPI2 <- inputDataPPI()
        incProgress(0.25, detail = "(Potential DMI)")
        Uni_DMI <- potentialDMIs()
        #PPI-DMI Mapping
        ########################################################################
        #names(PPI2) <- c("mProtein", "dProtein")
        incProgress(0.25, detail = "(Merge DMI & PPI)")
        predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
        Uni_predDMIs <- unique(predDMI)
        #names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
        predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
        adata$data$predDMI = predDMIs
        adata$data$predDMINR = unique(predDMIs[,c("mProtein","dProtein")])
        incProgress(0.5, detail = "(Complete)")
      })
    }
    predDMIs <- adata$data$predDMI
    print(head(predDMIs))
    return(predDMIs)
  })
  
  #Shows predicted DMIs in DataTable
  output$PredDMIs<-DT::renderDataTable({
    #Run only if Run button is active
    input$nrdmi
    if(input$run){
      style <- isolate(input$style)
      
      if(input$nrdmi){
        NRpredDMI <- predictedDMIs()
        print(head(NRpredDMI))
        NRpredDMI <- unique(NRpredDMI[,c("mProtein","dProtein")])
        print(head(NRpredDMI))
        formatStyle(datatable(NRpredDMI), columns = 1:2, color = "black")
        
      }else{
        formatStyle(datatable(predictedDMIs()), columns = 1:4, color = "black")
      }
    }
  })
  #*************************************************************************************
  
  #Step 3: Statistics
  ####################################################
  summaryStat <- eventReactive(input$run, {
    predDMI <- predictedDMIs()
    
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
      Datatype = factor(c("Motif","Domain","mProtein","dProtein")),
      Numbers = c(a,b,c,d)
    )
    
    p <- ggplot(data=statdmi, aes(x=Datatype, y=Numbers,fill=Datatype)) +
      geom_bar(colour="black", stat="identity") +
      guides(fill=FALSE)+
      scale_fill_manual(values = c("#9B59B6", "#85C1E9", "gold", "red"))
    
    p <- ggplotly(p)
    #p
    
    #barplot(x, main="Statistics of DMIs", col = colors)
    
    #legend("topright",
    #legend = c(paste("Motif=",a),paste("Domain=",b),paste("Motif containing Proteins=",c),paste("Domain containing Proteins=",d)), fill = colors)
    
  })
  output$plotbar <- renderPlotly({
    if(input$run){
      style <- isolate(input$style)
      
      summaryStat()
    }
  })
  #Step 4: Distribution of ELMs in the Predicted DMIs
  #####################################################
  #*************************************************************************************
  disELMs <- eventReactive(input$run, {
    #Read uploaded files
    GOterms <- read.csv("data/elm_goterms.tsv",header=TRUE,sep="\t")
    names(GOterms) <- c("ELM", "GO Term", "Biological Function")
    # Motif <- lowername(Motif)
    # Motif <- Motif[, c("accnum","motif")]
    # names(Motif) <- c("UniprotID","Motif")
    # Motif_NR<-unique(Motif)
    # #Rename the columns in two files
    # names(Motif_NR) <- c("mProtein", "Motif")
    # names(Domain) <- c("Motif", "Domain")
    # 
    # #Join/Merge two files based on Motif
    # Join <- merge( Motif_NR, Domain, by="Motif")
    # #print(Join)
    # names(Join) <- c("Motif", "Seq", "Domains")  #Change header of the output file
    # #Load mProtein_Motif_Domain file (result file from the previous code)
    # names(dProtein) <- c("Domains", "dProteins")
    # #joined both files based on Domains
    # DMI <- merge(Join, dProtein,by="Domains")
    # #Filtered unique DMIs
    # Uni_DMI <- unique(DMI)
    # #Named the header of output file
    # names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
    # #PPI-DMI Mapping
    # ########################################################################
    # names(PPI2) <- c("mProtein", "dProtein")
    # predDMI <- merge(PPI2, Uni_DMI, by= c("mProtein", "dProtein"))
    
    predDMIs <- predictedDMIs() # generatePredictedDMIs(input,adata$data)
    print(head(predDMIs))
    predDMI = predDMIs
    
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
    #!# This is not a pvalue! (No idea what it is)
    # pvalueelm <- round(df_pred2$Frequency/nrow(df_pred),2)
    # pvaluecol <- cbind(df_pred2,pvalueelm)
    propelm <- round(df_pred2$Frequency/sum(df_pred2$Frequency),2)
    propcol <- cbind(df_pred2,propelm)
    names(propcol) <- c("Motif", "Frequency", "Proportion")
    # if(input$DMIStrategy == c("elmcdom","elmcprot")){
    #   #!# Make and add GO toggle? Doesn't seem to work
    #   names(pvaluecol) <- c("ELM", "Frequency", "Pvalue")
    #   GeneOntology <- merge(pvaluecol,GOterms, by="ELM")
    #   return(GeneOntology)
    # }else{
    #   return(pvaluecol)
    # }
    return(propcol)
  })
  
  output$diselmsdata <-DT::renderDataTable({
    #Run only if Run button is active
    if(input$run){
      style <- isolate(input$style)
      disELMs()
      
    }
  })
  displotfunc <- function(){
    predDMI <- predictedDMIs()  # generatePredictedDMIs(input,adata$data)
    predDMIs = predDMI
    
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
      displotfunc()
    }
  })
  #Distribution of Domains in the Predicted DMIs
  #####################################################
  #*************************************************************************************
  disDom <- eventReactive(input$run, {
    
    predDMIs <- predictedDMIs()  # generatePredictedDMIs(input,adata$data)
    predDMI = predDMIs
    Uni_predDMIs <- unique(predDMI)
    names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
    #predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
    df_pred <- data.frame(Uni_predDMIs)["Domain"]
    names(df_pred) <- "Domain"
    print(df_pred)
    for (i in 1:length(df_pred)) {
      Matches <- count(df_pred)
      #names(Matches) <- "Frequency"
      #col <- cbind(df_pred,newcol)
      print(Matches)
      #
      
    }
    
    df_pred2 <- data.frame(Matches)[,(1:2)]
    names(df_pred2) <- c("Domain","Freq")
    print(df_pred2)
    Frequency <- df_pred2$Freq
    ELMs_names <- df_pred2$Domain
    df_pred2 <- data.frame(ELMs_names,Frequency)
    #pvalueelm <- round(df_pred2$Frequency/nrow(df_pred),2)
    #pvaluecol <- cbind(df_pred2,pvalueelm)
    propelm <- round(df_pred2$Frequency/sum(df_pred2$Frequency),2)
    propcol <- cbind(df_pred2,propelm)
    names(propcol) <- c("Domain", "Frequency", "Proportion")
    #names(pvaluecol) <- c("Domain", "Frequency", "Pvalue")
    propcol
    
  })
  
  output$disdomdata <-DT::renderDataTable({
    #Run only if Run button is active
    if(input$run){
      #Progress bar
      style <- isolate(input$style)
      disDom()
      
    }
  })
  disdomplotfunc <- function(){
    predDMIs <- predictedDMIs()  #generatePredictedDMIs(input,adata$data)
    predDMI = predDMIs
    
    Uni_predDMIs <- unique(predDMI)
    names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
    #predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
    df_pred <- data.frame(Uni_predDMIs)["Domain"]
    names(df_pred) <- "Domain"
    print(df_pred)
    for (i in 1:length(df_pred)) {
      Matches <- count(df_pred)
      #names(Matches) <- "Frequency"
      #col <- cbind(df_pred,newcol)
      print(Matches)
      #
      
    }
    df_pred2 <- data.frame(Matches)[,(1:2)]
    names(df_pred2) <- c("Domain","Freq")
    print(df_pred2)
    Frequency <- df_pred2$Freq
    dom_names <- df_pred2$Domain
    df_pred2 <- data.frame(dom_names,Frequency)
    disdomdmi <- data.frame(
      Datatype = factor(df_pred2$dom_names),
      Numbers = df_pred2$Frequency
    )
    
    p <- ggplot(data=disdomdmi, aes(x=Datatype, y=Numbers,fill=Datatype)) +
      geom_bar(colour="black", stat="identity") +
      guides(fill=FALSE)+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
    p <- ggplotly(p)
    
    
  }
  
  output$disdomchart <- renderPlotly({
    #Run only if Run button is active
    if(input$godisd){
      #Progress bar
      style <- isolate(input$style)
      disdomplotfunc()
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
    Domain <- inputDataMotifDomain()
    dProtein <- inputDatadomain()
    Motif <- inputDataMotif()
    PPI <- inputDataPPI()
    
    Motif_NR<-unique(Motif)
    
    # #Join/Merge two files based on Motif
    Join <- merge( Motif_NR, Domain, by="Motif")
    #print(Join)
    # names(Join) <- c("Motif", "Seq", "Domain")  #Change header of the output file
    # #Load mProtein_Motif_Domain file (result file from the previous code)
    # names(dProtein) <- c("Domains", "dProteins")
    # #joined both files based on Domains
    DMI <- merge(Join, dProtein,by="Domain")
    # #Filtered unique DMIs
    Uni_DMI <- unique(DMI)
    # #Named the header of output file
    # names(PPI) <- c("mProtein", "dProtein")
    # names(Uni_DMI) <- c("Domain", "Motif", "mProtein", "dProtein")
    predDMI <- merge(PPI, Uni_DMI, by= c("mProtein", "dProtein"))
    Uni_predDMIs <- unique(predDMI)
    # names(Uni_predDMIs) <- c("mProtein", "dProtein", "Domain", "Motif")
    predDMIs <- Uni_predDMIs[, c("mProtein","Motif", "Domain", "dProtein")]
    
    # Reduce Uni_DMI and predDMIs to unique mProtein-dProtein pairs    
    Uni_DMI <- unique(DMI[,c("mProtein","dProtein")])
    predDMIs <- unique(predDMI[, c("mProtein","dProtein")])
    
    
    PPI <- inputDataPPI()
    PPI_data <- unique(PPI)
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
    showNotification("Performing the randomisations", type = "message", duration = 5)
    data <- list()
    withProgress(message = 'Performing randomisations', detail = 0, value = 0, {
      for (j in 1:input$shufflenum) {
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
        incProgress(1/input$shufflenum, detail = j)
      }
    })
    #rPPI-DMI Mapping                                                               
    #################################################################################
    showNotification(paste(input$shufflenum,"random PPI datasets have been created"), type = "message", duration = 5)
    showNotification("Now predicting DMIs from the random PPI data", type = "message", closeButton = TRUE,duration = 15)
    m <- data.frame()
    withProgress(message = 'Predicting random DMI', detail = 0, value = 0, {
      for (i in 1:input$shufflenum) {
        rPPI <- data[[i]]
        names(rPPI)<-c("mProtein", "dProtein")
        DMI_rPPI <- merge(Uni_DMI, rPPI, by= c("mProtein", "dProtein"))
        DMI_rPPI <- unique(DMI_rPPI[,c("mProtein", "dProtein")])
        Matches <- nrow(DMI_rPPI)
        if(devmode){ print(Matches) }
        dmatch <- data.frame(Matches)
        m=rbind(m,dmatch)
        row.names(m) <- NULL
        incProgress(1/input$shufflenum, detail = i)
      }
    })
    m
  })
  #*************************************************************************************
  
  #Step 6: Histogram
  ####################################################
  
  output$histogram <- renderPlot({
    
    if(input$run){
      style <- isolate(input$style)
      showNotification("Creating plot - this may take a while! Nothing will respond while it's being calculated (e.g. network tab)", type = "warning", duration = NULL, id="plotwarning")    
      #withProgress(message = 'Creating Plot', detail = "This may take a while!!. Nothing will respond while it's being calculated (e.g. network tab)", value = 0, {
      rPPIDMI()
      plotInput()
      removeNotification("plotwarning")
    }
  })
  #Function to generate Histogram
  plotInput <- function(){
    #File upload check
    Domain <- inputDataMotifDomain()
    dProtein <- inputDatadomain()
    Motif <- inputDataMotif()
    PPI <- inputDataPPI()
    
    #dirName <- paste0("RandomFiles_", strsplit(as.character(PPIFile$name), '.csv'))
    x <- rPPIDMI()
    predictedDMImutlilink = predictedDMIs()
    predDMIs = unique(predictedDMImutlilink[,c("mProtein","dProtein")])
    predDMInum = nrow(predDMIs)
    
    names(x) <- "values"
    #bins<- seq(min(x), max(x), length.out = input$bins + 1)
    par(bg=input$col2)
    xlimmax = max(c(input$xlimend,max(x)*1.1,predDMInum*1.1))
    bins<- seq(0, xlimmax+input$binwidth, input$binwidth)
    # if(max(bins) < xlimmax){
    #   bins = c(bins,max(bins)+input$binwidth)
    # }
    #ymax = max(hist(x$values, breaks=bins, plot=FALSE)$counts) * 1.5   
    h<- hist(x$values, breaks=bins, col = input$col, border = 'black', main=input$text3, ylab=input$text2,
             xlab=input$text, xlim = c( input$xlimstart,  xlimmax), labels = input$barlabel, cex.main=1.5, cex.lab=1.5,cex.axis=1.5)

    xfit <- seq(min(x$values),max(x$values),length=40)
    yfit <- dnorm(xfit, mean=mean(x$values),sd=sd(x$values))
    yfit<- yfit*diff(h$mids[1:2])*length(x$values)
    lines(xfit,yfit)

    ob_fdr <- x[x$values <= predDMInum, ]
    #axis(side=3, lwd = 0, lwd.ticks = 4, at=nrow(predDMIs), lend=1, labels = FALSE, tcl=5, font=2, col = "black", padj = 0, lty = 3)
    #shows the observed value
    mtext(paste("Observed value: ", predDMInum), side = 3, at=predDMInum, font = 4)
    pvalue = length(x[x >= predDMInum])/input$shufflenum
    
    pplace = xlimmax
    if(predDMInum > xlimmax/2){ pplace = predDMInum/2 }
    
    if(pvalue == 0){ 
      mtext(paste0("P-value is: < ", 1/input$shufflenum), side = 3, at=pplace, font = 4, col = "red")
      pvalue <-  paste0("<b>P-value: </b> < ", 1/input$shufflenum)
    }else{
      #mtext(paste0("P-value is: ", pvalue), side = 3, at=predDMInum+90, font = 4, col = "red")
      mtext(paste0("P-value is: ", pvalue), side = 3, at=pplace, font = 4, col = "red")
      pvalue <-  paste0("<b>P-value: </b>", pvalue)
    }
    #mtext(paste0("Mean is: ", mean(x$values)), side = 3, at=mean(x$values), font = 4, col="red")
    #mtext(paste0("FDR is: ", mean(x$values)/nrow(predDMIs)), side = 3, at=50, font = 4, col= "red")
    #points arrow on the observed value
    arrows(predDMInum, 480, predDMInum, 0, lwd = 2, col = "black", length = 0.1, lty = 3)
    obsval <- paste0("<b>Observed NR DMI: </b>", predDMInum)
    obsdmi <- paste0("<b>Observed DMI (redundant): </b>", nrow(predictedDMImutlilink))
    meanvalue <- paste0("<b>Mean random DMI: </b>", round(mean(x$values),2))
    Escore <- paste0("<b>Enrichment score (E-score): </b>", round(predDMInum/mean(x$values),2))
    FDR <- paste0("<b>False Discovery Rate (FDR): </b>", round(mean(ob_fdr)/predDMInum,2))
    output$summary <- renderUI({
      
      HTML(paste("<font color=\"#FF0000\"><b>Summary of Histogram</b></font>", pvalue, obsdmi, obsval, meanvalue, Escore, FDR, sep = '<hr/>'))
      
    })
  }
  
  #?# What does this function actually do?!
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
  
  #Step 7: DMI Network
  ####################################################
  mynetwork <- function(){
    set.seed(20)
    if(input$run){
      style <- isolate(input$style)
      
      predDMIs <- predictedDMIs() # generatePredictedDMIs(input,adata$data)
      predDMI = predDMIs
      
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
        
        #   }
        # }
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
      png(file, width = input$width, height = input$height, units = "px", pointsize = 12, res = 300)
      plotInput()
      dev.off()
    }
  )
  #download button to download the plot
  output$downloadrandom <- downloadHandler(
    filename = "predDMIs.csv",
    content = function(file) {
      write.csv(predictedDMIs(), file)
    }
  )
})
##############################
# End of Server code.
##############################

#*********************************************************************************************************
#*********************************************************************************************************

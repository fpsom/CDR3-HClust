library(shinyFiles)
library("shinyBS")
library("DT")
library(dplyr)
library(plotly)
library(xtable)


source("helpers.R")



shinyServer( 
  function(session, input, output) { 
    
    ######################################## initialize global variables  ########################################
    # A list with the similarity groups
    sim = list("F","W",c("A","I","L","V"),c("M","C"),"P","G","Y",c("T","S"),c("H","K","R"),c("E","D"),c("Q","N"))
    # Naming the group of similarities
    names(sim) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
    # A table with the letters
    let = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y") # The letters matrix
    # An empty vector wich will store the level for every cluster
    clep = vector('numeric')
    br = 0  # Initial value of branch
    cl = 0  # Initial value of new clusters
    met = 0 # Initial value of level capacity counter
    ep = 1  # Initial value of level
    d = 0
    nn = FALSE # Initial value for the condition sumper < endper%
    endper = 90 # Set the percentage, which ends the programm 
    listxx = list() # Initialize a list for saving the permat of all branches
    listyy = list()
    dfsum = data.frame(sumper = numeric(0),sumper2 = numeric(0),branch = numeric(0), len = numeric(0))
    dfsd = data.frame(Average_Identity_Value = numeric(0),Identity_Standar_Deviation = numeric(0),Average_Similarity_Value = numeric(0), Similarity_Standar_Deviation = numeric(0)) 
    ############################################# Select Datasets #############################################
    
    # dir
    shinyDirChoose(input, 'dir', roots = c(home = '.'), filetypes = c('', 'csv'))
    
    dir <- reactive(input$dir)
    output$dir <- renderPrint(dir())
    
    # path
    path <- reactive({
      home <- normalizePath(".")
      file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
    })
    
    # files
    output$uiDatasets <- renderUI({
      if (is.null(dir())) return()
      checkboxGroupInput(inputId = "Dataset", label = "Select Datasets", inline=TRUE, choices = list.files(path()))
    })
    
    observeEvent(is.null(input$begin),{
      if (is.null(input$Dataset)) return()
      udata=read.csv(paste0("data/",input$Dataset), header = TRUE, sep = ";")
      udata$AA.JUNCTION = as.character(udata$AA.JUNCTION)
      udata$clusters = 0 # Initialiaze the column clusters with 0
      udata$level.0 = 0 # Initialize the column of cl.0 with 0
      udata$temp= 0 # Creating a temp column with 0
      mm = 0
      listq = list("dfsum" = dfsum,"list" = listxx, "listn" = listyy,"udata" = udata,"permat"= NA, "persim" = NA, "br" = br, "cl" = cl, "met" = met, "ep" = ep , "clep" = clep, "nn" = nn, "sumper" = NA,"sumper2" = NA, "ela" = NA, "cel" = NA,"endper" = endper)
      lastlist = Matrices(listq,FALSE,let,sim,d)
      # The final name of udata data frame 
      df = lastlist$udata
      # A list with the permat matrix for every branch
      perlist = lastlist$list
      persimlist = lastlist$listn
      # Create a dataframe with all clusters and their identity and similarity percentage 
      ff = lastlist$dfsum
      ff$level[1]= 0
      ff$level[2:(length(lastlist$clep[ff$branch]) +1 )] = lastlist$clep[ff$branch]  #without leaves
    
      Clus = as.data.frame(matrix(100, ncol = 3, nrow = max(df$clusters)+1))
      names(Clus) = c("ClusterId","Identity","Similarity")
      Clus$ClusterId = 0:max(df$clusters)
      Clus$seqnum = 1
      Clus$level = c(0,lastlist$clep)
      for(i in 1:length(ff$branch) ){
        ll = which(Clus$ClusterId == ff$branch[i])
        Clus$Identity[ll] = ff$sumper[i]
        Clus$Similarity[ll] = ff$sumper2[i]
        Clus$seqnum[ll] = ff$len[i]
      }
    
      mm = 1
    
      output$message <- renderText({
        if(is.null(mm)) return()
        "Data available from now on"
      })
    
      output$tre <- renderGrViz({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        Den(input$tree_level,df,lastlist)
      })
    
      output$coltre <- renderCollapsibleTree({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        cc(df,Clus)
      })
    
      output$logolev <- renderPlot({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        plot(LogoLev(input$logo_level,lastlist,df))
      })
    
      output$downloadLogoLevel <- downloadHandler(
        filename = function(){paste0("logo_level_",input$logo_level,".png")},
        content = function(file) {
          png(file,width=1000, height=550)
          plot(LogoLev(input$logo_level,lastlist,df))
          dev.off()
        }) 
    
      output$logocl <- renderPlot({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        plot(LogoCl(input$logo_cluster,lastlist,df))
      }) 
    
      output$downloadLogoCluster <- downloadHandler( 
        filename = function(){paste0("logo_cluster_",input$logo_cluster,".png")},
        content = function(file) {
          png(file,width=1000, height=550)
          plot(LogoCl(input$logo_cluster,lastlist,df))
          dev.off()
      }) 
    
      output$barlev <- renderPlot({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        BarLev(input$bar_level,lastlist,perlist,persimlist,Clus,let,sim,input$barl_style)
      })
    
      output$downloadBarLevel <- downloadHandler(
        filename = function(){paste0("barPlot_level_",input$bar_level,".png")},
        content = function(file) {
          png(file,width=1000, height=550)
          BarLev(input$bar_level,lastlist,perlist,persimlist,Clus,let,sim,input$barl_style)
          dev.off()
      }) 
    
      output$barcl <- renderPlot({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        BarCl(input$bar_cluster,perlist,persimlist,Clus,let,sim,input$barcl_style)
      })
    
      output$downloadBarCluster <- downloadHandler(
        filename = function(){paste0("barPlot_cluster_",input$bar_cluster,".png")},
        content = function(file) {
            png(file,width=1000, height=550)
            BarCl(input$bar_cluster,perlist,persimlist,Clus,let,sim,input$barcl_style)
            dev.off()
      }) 
    
      output$aminolev <- renderTable({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        AminoLev(input$amino_level,lastlist,df,Clus)
      })
    
      output$aminocl <- renderTable({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        AminoCl(input$amino_cluster,lastlist,df)
      })
    
      output$pin <- renderTable({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        EmPin(lastlist,Clus,dfsd)
      })
    
      output$downloadallEmPin <- downloadHandler(
        filename = function(){"identity_per_level_table.txt"},
        content = function(file) {
          write.table(EmPin(lastlist,Clus,dfsd),file,sep="\t")
      }) 
    
      output$idenaminolev <- renderTable({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        idenlev(input$aminoide_level,Clus,input$aminoide_lever)
      })
    
      output$idenaminocl <- renderTable({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        idencl(input$aminoide_cluster,Clus,input$aminoide_clver)
      })
    
      output$opt1lev <- renderText({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        Opt(AminoLev(input$opt1_level,lastlist,df,Clus))
      })
    
      output$opt1cl <- renderText({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        Opt(AminoCl(input$opt1_cluster,lastlist,df))
      })
    
      output$opt2lev <- renderText({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        Opt2(AminoLev(input$opt2_level,lastlist,df,Clus),sim)
      })
    
      output$opt2cl <- renderText({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        Opt2(AminoCl(input$opt2_cluster,lastlist,df),sim)
      })
    
      output$idplot <- renderPlot({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        plot(Id(ff,input$ident))
      })
    
      output$downloadIdentityPlot <- downloadHandler(
        filename = function(){paste0(input$ident,"_plot",".png")},
        content = function(file) {
          png(file,width=1000, height=550)
          plot(Id(ff,input$ident))
          dev.off()
      })
      
      output$net <- renderPlot({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        Netw(input$net_level,input$thr_value,input$thr_style,input$net_style,df,lastlist,Clus)
      })
      
      output$netInt <- renderNdtvAnimationWidget({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        
        NetwInt(input$netInt_level,input$thrInt_value,input$thrInt_style,input$netInt_style,df,lastlist,Clus)
        
        #render.d3movie(nn, usearrows = F, displaylabels = TRUE, bg="#ffffff", 
        #               vertex.border="#ffffff", vertex.col =  nn %v% "color",
        #               vertex.cex = nn %v% "size", label = nn %v% "label", edge.lwd = nn %e% "width" , edge.col = nn %e% "color",output.mode = 'htmlWidget',script.type= 'embedded')
      })
    
   })    
 })

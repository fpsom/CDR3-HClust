source("helpers.R")



shinyServer( 
  function(session, input, output) { 
    
    ######################################## initialize global variables  ########################################
    # A list with the similarity groups
    sim = list("F","W",c("A","I","L","V"),c("M","C"),"P","G","Y",c("T","S"),c("H","K","R"),c("E","D"),c("Q","N"))
    # Naming the group of similarities
    names(sim) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
    altsim = list("F","W",c("A","I","L","V"),c("M","C"),"P","G","Y",c("T","S"),c("H","K","R"),c("E","D"),c("Q","N"))
    # Naming the group of similarities
    names(altsim) = c("f","w","a","s","p","g","y","h","b","c","m")
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
    dfsd = data.frame(ff = character(0),Average_Identity_Value = numeric(0),Identity_Standar_Deviation = numeric(0),Average_Similarity_Value = numeric(0), Similarity_Standar_Deviation = numeric(0), stringsAsFactors = FALSE) 
    ggdf = data.frame(branch = numeric(0), len = numeric(0))
    last = 0
    flagtic = TRUE
    progend = FALSE
    listax = list()
     ############################################# Select Datasets #############################################
    
    dfup <- reactive({
      inFiles <- input$fileIn
      print(inFiles)
      dfup <- data.frame()
      if (is.null(inFiles))
        return(NULL)
      for (i in seq_along(inFiles$datapath)) {
        tmp <- read.csv(inFiles$datapath[i], header = TRUE, sep = ";")  
        dfup <- rbind(dfup, tmp)
      }
      dfup
      
    })
    output$tbl <- DT::renderDataTable(
      if(input$fdet == TRUE){
        dfup()
      }
    )
    output$tbl2 <- DT::renderDataTable(
      if(input$fdet == TRUE){
        input$fileIn
      }
    )
    output$results = renderPrint({
    })
    
    observeEvent(is.null(input$begin),{
      if (is.null(input$fileIn)) return()
      hide("begin")
      udata = dfup()
      udata$AA.JUNCTION = as.character(udata$AA.JUNCTION)
      udata$clusters = 0 # Initialiaze the column clusters with 0
      udata$level.0 = 0 # Initialize the column of cl.0 with 0
      udata$temp= 0 # Creating a temp column with 0
      mm = 0
      listax$temp = 1:length(udata$AA.JUNCTION)
      names(listax)[length(listax)] = sprintf('cl.%d', 0) 
      listq = list("listax" = listax,"ggdf" = ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy,"udata" = udata,"permat"= NA, "persim" = NA, "br" = br, "cl" = cl, "met" = met, "ep" = ep , "clep" = clep, "nn" = nn, "sumper" = NA,"sumper2" = NA, "ela" = NA, "cel" = NA,"endper" = endper, "last" = last,"progend" = progend, "leaf" = FALSE)
      listb = listq
      
      #logfile
      logFile = paste0(getwd(),"/log_file ",trunc(as.numeric(Sys.time())),".txt")
      cat(paste0("Function","\t","Num of input rows","\t","Num of input columns","\t","Duration"), file=logFile, append=FALSE, sep = "\n")
      
      tic()
      while (progend == FALSE){
        lista = Matrices(listb,leaf,let,sim,d,input$algorithm_type,input$column_num,logFile)
        progend = lista$progend
        if(progend == TRUE){
          break
        }
        listb = Choice(lista,leaf,let,sim,d,input$algorithm_type,input$column_num,logFile)
        progend = listb$progend
      }
      en = toc(quiet = TRUE)
      cat(paste0("lastlist","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
      lastlist = listb
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
      Clus$seqnum[1] = nrow(df) 
      Clus$seqnum[2:(max(df$clusters)+1)] = lastlist$ggdf$len
      Clus$level = c(0,lastlist$clep)
      for(i in 1:length(ff$branch) ){
        ll = which(Clus$ClusterId == ff$branch[i])
        Clus$Identity[ll] = ff$sumper[i]
        Clus$Similarity[ll] = ff$sumper2[i]
      }
      
      
      enthr = input$thr_en
      if(enthr == TRUE){
        lev = max(na.omit(lastlist$clep))
        threshold1 = input$thr_1
        threshold2 = input$thr_2
        algo = input$algorithm_type
        df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
        df_args <- c(df[str_which(names(df), "level.")], sep="/")
        if(lev == max(na.omit(lastlist$clep))){
          df$pathString<- do.call(paste, df_args)
          kk = df$pathString
          for(i in 1:length(kk)){
            temp = str_locate(kk[i],"/NA")
            if(is.na(temp[1]) == FALSE){
              temp2 = str_sub(kk[i], 1, temp[1]-1);
              kk[i] = temp2 
            }
          } 
        }else{
          df$pathString<- do.call(paste, df_args)
          kk = df$pathString
          gg =as.data.frame(str_locate_all(kk,"/"))
          tem = 1
          for (i in 1:length(kk)){
            kk[i] = str_sub(kk[i],1,gg[,tem][lev+1]-1)
            tem = tem + 2
          }
          for(i in 1:length(kk)){
            temp = str_locate(kk[i],"/NA")
            if(is.na(temp[1]) == FALSE){
              temp2 = str_sub(kk[i], 1, temp[1]-1);
              kk[i] = temp2 
            }
          } 
        }
        df$pathString = kk
        x <- ToDataFrameTree(df, "pathstring")
        xN <- as.Node(x)
        
        matches <- regmatches(unique(x$pathString), gregexpr("[[:digit:]]+", unique(x$pathString)))
        arxpath = unique(x$pathString)
        for(i in 1:length(matches)){
          ggg = as.numeric(unlist(matches[i])) +1
          for(j in 1:(length(ggg)-1)){
            if(algo == "Identity"){
              ddd1 = ((Clus$seqnum[ggg[j]] - Clus$seqnum[ggg[j+1]]) / Clus$seqnum[ggg[j]]) * 100
              ddd2 =(Clus$Identity[ggg[j+1]] - Clus$Identity[ggg[j]]) * 0.5 + (Clus$Similarity[ggg[j+1]] - Clus$Similarity[ggg[j]])* 0.5
              #ddd1 > threshold1 && ddd2 > threshold2
              if(ddd1 > threshold1 && ddd2 > threshold2){
                deik = which(x$pathString == arxpath[i])
                tem = str_locate(arxpath[i],as.character(ggg[j]-1))
                fftemp = str_sub(x$pathString[deik[1]],1,tem[2])
                ff2 = str_sub(x$pathString[deik[1]],tem[1],tem[2])
                ff3 = as.numeric(ff2)
                x$clusters[deik] = ff3
                ff4 = lastlist$clep[ff3]
                ff5 = str_which(names(x),sprintf("level.%d",ff4))
                ff6 = str_which(names(x), "level.")
                if(is.na(ff5[1]) == FALSE){
                  x[deik,(ff5[1]+1):max(ff6)] = NA
                }
                x$pathString[deik] = fftemp
                break()
              }
            }else if(algo == "Similarity"){
              ddd1 = ((Clus$seqnum[ggg[j]] - Clus$seqnum[ggg[j+1]]) / Clus$seqnum[ggg[j]]) * 100
              ddd2 = Clus$Similarity[ggg[j+1]] - Clus$Similarity[ggg[j]]
              if(ddd1 > threshold1 && ddd2 > threshold2){
                deik = which(x$pathString == arxpath[i])
                tem = str_locate(arxpath[i],as.character(ggg[j]-1))
                fftemp = str_sub(x$pathString[deik[1]],1,tem[2])
                ff2 = str_sub(x$pathString[deik[1]],tem[1],tem[2])
                ff3 = as.numeric(ff2)
                x$clusters[deik] = ff3
                ff4 = lastlist$clep[ff3]
                ff5 = str_which(names(x),sprintf("level.%d",ff4))
                ff6 = str_which(names(x), "level.")
                if(is.na(ff5[1]) == FALSE){
                  x[deik,(ff5[1]+1):max(ff6)] = NA
                }
                x$pathString[deik] = fftemp
                break()
              }
            }
          }
        }
        xN <- as.Node(x)
        bo = which(lastlist$clep == lev)
        levtel = lastlist$clep
        j = 1
        jfjf = vector()
        for(i in 1:max(bo)){
          temp1 = as.numeric(names(FindNode(xN,(sprintf("%d",i)))$children))
          if(length(temp1) == 1){
            if(temp1[1] %% 2 == 1){
              ff7 = temp1[1] + 1
            }else{
              ff7 = temp1[1] - 1
            }
            ff8 = lastlist$clep[ff7]
            ff9 = str_which(names(x),sprintf("level.%d",ff8))
            deik2 = which(df[ff9] == ff7)
            x[deik2,ff9] = df[deik2,ff9]
            tem = str_locate(df$pathString[deik2[1]],as.character(ff7))
            fftemp = str_sub(df$pathString[deik2[1]],1,tem[2])
            x$pathString[deik2] = fftemp
            jfjf[j] = ff7
            j = j+1
          }
          if(is.null(FindNode(xN,(sprintf("%d",i)))$isLeaf) && is.element(i, jfjf) == FALSE){
            levtel[i] = NA
            Clus[i+1,] = NA
            ff[i+1,] = NA
          }
        }
        xN <- as.Node(x)
        rm(lev)
        lastlist$udata = x
        df = x
        lastlist$clep = levtel
      }
      
      mm = 1
      #fname = input$fileIn$name
      
      output$message <- renderText({
        if(is.null(mm)) return()
        "Data available from now on"
      })
      
      
      #val2<-lastlist
      
      #makeReactiveBinding("val2")
      
      # Save extra values in state$values when we bookmark
      onBookmark(function(state) {
        print("dfdfadsfffffffffffffffffffffffffff")
        state$values$val2<-lastlist
        print(state$values$val2)
      })
      
      # Read values from state$values when we restore
      onRestore(function(state) {
        val2<<-state$values$val2
        print("###############################")
        print(val2)
      })
 
      
      observeEvent(input$showPreviousSessions, {
        
        output$uiPreviousSessions <- renderUI({
          if (is.null(list.files('shiny_bookmarks'))) return()
          dates_of_files=lapply(list.files('shiny_bookmarks'),function(x) file.mtime(paste0("shiny_bookmarks/",x)))
          
          print(dates_of_files)
          checkboxGroupInput(inputId = "fileIn", label = "Select Datasets", inline=TRUE, choices = list.files('shiny_bookmarks'))
        })
        
        output$uiLoad <- renderUI({
          if (is.null(input$fileIn)) return()
          link=paste0("http://127.0.0.1:5211/?_state_id_=",input$fileIn)
          
          
          #actionButton("load", "Load",onclick="window.open('http://127.0.0.1:7853/?_state_id_=e96cfb203652fd1a', '_blank')")
          
          helpText(a("Click Here to load session", href=link, target="_blank")
                   
          )
        })
      })
      
    
      
      output$tre <- renderGrViz({
        if (is.null(input$fileIn)) return()
        Den(input$tree_level,df,lastlist,Clus,flagtic,input$algorithm_type,threshold1,threshold2,input$thr_en,logFile)
      })
      
      output$coltre <- renderCollapsibleTree({
        if (is.null(input$fileIn)) return()
        cc(df,Clus,flagtic,logFile)
      })
    
      output$logolev <- renderPlot({
        if (is.null(input$fileIn)) return()
        lole = LogoLev(input$logo_level,lastlist,df,flagtic,logFile)
        if (is.null(lole)){
          plot(1, type="n", axes=F, xlab="", ylab="")
        }else{
          plot(lole)
        }
      })
    
      output$downloadLogoLevel <- downloadHandler(
        filename = function(){paste0("logo_level_",input$logo_level,".png")},
        content = function(file) {
          png(file,width=1000, height=550)
          plot(LogoLev(input$logo_level,lastlist,df,flagtic,logFile))
          dev.off()
        }) 
    
      output$logocl <- renderPlot({
        if (is.null(input$fileIn)) return()
        locl = LogoCl(input$logo_cluster,lastlist,df,flagtic,logFile)
        if (is.null(locl)){
          plot(1, type="n", axes=F, xlab="", ylab="")
        }else{
          plot(locl)
        }
      }) 
    
      output$downloadLogoCluster <- downloadHandler( 
        filename = function(){paste0("logo_cluster_",input$logo_cluster,".png")},
        content = function(file) {
          png(file,width=1000, height=550)
          plot(LogoCl(input$logo_cluster,lastlist,df,flagtic,logFile))
          dev.off()
      }) 
    
      output$barlev <- renderPlot({
        if (is.null(input$fileIn)) return()
        BarLev(input$bar_level,lastlist,perlist,persimlist,Clus,let,sim,input$barl_style,flagtic,logFile)
      })
    
      output$downloadBarLevel <- downloadHandler(
        filename = function(){paste0("barPlot_level_",input$bar_level,".png")},
        content = function(file) {
          png(file,width=1000, height=550)
          BarLev(input$bar_level,lastlist,perlist,persimlist,Clus,let,sim,input$barl_style,flagtic,logFile)
          dev.off()
      }) 
    
      output$barcl <- renderPlot({
        if (is.null(input$fileIn)) return()
        BarCl(input$bar_cluster,perlist,persimlist,Clus,let,sim,input$barcl_style,lastlist,flagtic,logFile)
      })
    
      output$downloadBarCluster <- downloadHandler(
        filename = function(){paste0("barPlot_cluster_",input$bar_cluster,".png")},
        content = function(file) {
            png(file,width=1000, height=550)
            BarCl(input$bar_cluster,perlist,persimlist,Clus,let,sim,input$barcl_style,lastlist,flagtic,logFile)
            dev.off()
      }) 
    
      output$aminolev <- renderTable({
        if (is.null(input$fileIn)) return()
        AminoLev(input$amino_level,lastlist,df,Clus,flagtic,logFile)
      })
    
      output$aminocl <- renderTable({
        if (is.null(input$fileIn)) return()
        AminoCl(input$amino_cluster,lastlist,df,flagtic,logFile)
      })
    
      output$pin <- renderTable({
        if (is.null(input$fileIn)) return()
        EmPin(lastlist,Clus,dfsd,flagtic,logFile)
      })
    
      output$downloadallEmPin <- downloadHandler(
        filename = function(){"identity_per_level_table.txt"},
        content = function(file) {
          write.table(EmPin(lastlist,Clus,dfsd,flagtic,logFile),file,sep="\t")
      }) 
    
      output$idenaminolev <- renderTable({
        if (is.null(input$fileIn)) return()
        idenlev(input$aminoide_level,Clus,input$aminoide_lever,flagtic,logFile)
      })
    
      output$idenaminocl <- renderTable({
        if (is.null(input$fileIn)) return()
        idencl(input$aminoide_cluster,Clus,input$aminoide_clver,flagtic,logFile)
      })
    
      output$opt1lev <- renderText({
        if (is.null(input$fileIn)) return()
        Opt(AminoLev(input$opt_level,lastlist,df,Clus,flagtic,logFile),flagtic,logFile)
      })
      
      output$opt1Newlev <- renderText({
        if (is.null(input$fileIn)) return()
        OptNew(AminoLev(input$opt_level,lastlist,df,Clus,flagtic,logFile),logFile)
      })
      
      output$opt1cl <- renderText({
        if (is.null(input$fileIn)) return()
        Opt(AminoCl(input$opt_cluster,lastlist,df,flagtic,logFile),flagtic,logFile)
      })
      
      output$opt1Newcl <- renderText({
        if (is.null(input$fileIn)) return()
        OptNew(AminoCl(input$opt_cluster,lastlist,df,flagtic,logFile),logFile)
      })
    
      output$opt2lev <- renderText({
        if (is.null(input$fileIn)) return()
        Opt2(AminoLev(input$opt_level,lastlist,df,Clus,flagtic,logFile),sim,FALSE,altsim,flagtic,logFile)
      })
      
      output$opt2Newlev <- renderText({
        if (is.null(input$fileIn)) return()
        Opt2New(AminoLev(input$opt_level,lastlist,df,Clus,flagtic,logFile),sim,logFile)
      })
    
      output$opt2cl <- renderText({
        if (is.null(input$fileIn)) return()
        Opt2(AminoCl(input$opt_cluster,lastlist,df,flagtic,logFile),sim,FALSE,altsim,flagtic,logFile)
      })
      
      output$opt2Newcl <- renderText({
        if (is.null(input$fileIn)) return()
        Opt2New(AminoCl(input$opt_cluster,lastlist,df,flagtic,logFile),sim,logFile)
      })
    
      output$idplot <- renderPlot({
        if (is.null(input$fileIn)) return()
        plot(Id(ff,input$ident,flagtic,logFile))
      })
    
      output$downloadIdentityPlot <- downloadHandler(
        filename = function(){paste0(input$ident,"_plot",".png")},
        content = function(file) {
          png(file,width=1000, height=550)
          plot(Id(ff,input$ident,flagtic,logFile))
          dev.off()
      })
      
      me <- reactiveVal(1)

      observeEvent(input$Apply, {
        netlist = Netw(input$net_level,input$thr_value,input$thr_style,input$net_style,df,lastlist,Clus,sim,altsim,input$table_show,FALSE,input$sil,flagtic,input$algorithm_type,input$thr_1,input$thr_2,input$thr_en,logFile)
        pal1 <- rainbow(6, alpha=1) 
        netyp = input$net_style
        net0.copy = netlist$net0.copy
        net1.copy = netlist$net1.copy
        thrtyp = netlist$thrtyp
        thrtyp2 = netlist$thrtyp2
        te = max(na.omit(as.vector(thrtyp)))
        bb = netlist$bb
        colrs = netlist$colrs
        output$net <- renderPlot({
          if (is.null(input$fileIn)) return()
          if(netyp == "whole"){
            par(mfrow = c(1,2))
            plot(net0.copy, edge.color=na.omit(pal1[( E(net0.copy)$width %/% 1) +1]), edge.curved=.1, vertex.label.color = "black") #plot the network graph
            legend("topleft", inset=c(0.1,0.2), colnames(df[str_which(names(df), "level.")])[1:(max(bb)+1)], pch=21,
                   col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
            legend("topright", inset=c(0.1,0.2), c("0-1","1-2","2-3","3-4","4-5","5"), pch=21,
                   col="#777777", pt.bg=pal1, pt.cex=2, cex=.8, bty="n", ncol=1,title = "Relationship strength")
            matrix.heatmap(thrtyp)
          }else{
            par(mfrow = c(1,2))
            plot(net1.copy, edge.color=na.omit(pal1[( E(net1.copy)$width %/% 1) +1]), edge.curved=.1, vertex.label.color = "black") #plot the network graph
            legend("topleft", inset=c(0.1,0.2), colnames(df[str_which(names(df), "level.")])[1:(max(bb)+1)], pch=21,
                   col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
            legend("topright", inset=c(0.1,0.2), c("0-1","1-2","2-3","3-4","4-5","5"), pch=21,
                   col="#777777", pt.bg=pal1, pt.cex=2, cex=.8, bty="n", ncol=1,title = "Relationship strength")
            matrix.heatmap(thrtyp2)
          }
        })
        
        output$dis_max <- renderText({
          me(te)
          if (input$thr_style == "Distance"){
            paste("The maximum distance is",te)
          }else if (input$thr_style == "StrIdentSimilarity"){
            paste("The maximum Identity Similarity is",te )
          }else if (input$thr_style == "StrGroupSimilarity"){
            paste("The maximum Group Similarity is",te)
          }else if (input$thr_style == "FinalDis"){
            paste("The maximum final distance is",te)
          }
        })
        
        output$net_table <- renderDataTable({
          if(input$table_show == TRUE){
            if(netyp == "whole"){
              thrtyp
            }else{
              thrtyp2
            }
          }
        })
        reset("Apply")
      })
      
      
      
      inline = function (x) {
        tags$div(style="display:inline-block;", x)
      }
      
        
    output$g1 <- renderUI({
      conditionalPanel(
        condition = "input.Tree % 2 == 1",
        br(),
        br(),
        numericInput("tree_level", "Select max level of the tree:", 5,  min = 1, max = max(na.omit(lastlist$clep)), width="140px"),
        br(),
        grVizOutput("tre", width = "1500px", height = "1000px")
      )
    })
    
    output$g2 <- renderUI({
      conditionalPanel(
      condition = "input.Coltree % 2 == 1",
      br(),
      br(),
      collapsibleTreeOutput("coltre", width = "1500px", height = "1000px")
    )})
    
    output$g3 <- renderUI({
      conditionalPanel(
        condition = "input.Logo % 2 == 1",
        br(),
        br(),
        numericInput("logo_level", "Select level to show:", 5,  min = 1, max = max(na.omit(lastlist$clep)), width="140px"),
        br(),
        plotOutput("logolev"),
        br(),
        downloadButton("downloadLogoLevel")  
      )
    })
    
    output$g4 <- renderUI({
      conditionalPanel(
        condition = "input.Logo % 2 == 1",
        br(),
        br(),
        numericInput("logo_cluster", "Select cluster to show:", 5,  min = 1, max = max(df$clusters), width="140px"),
        br(),
        plotOutput("logocl"),
        br(),
        downloadButton("downloadLogoCluster")
      )
    })
    
    output$g5 <- renderUI({
      conditionalPanel(
        condition = "input.Barplot % 2 == 1",
        br(),
        br(),
        numericInput("bar_level", "Select level to show:", 5,  min = 1, max = max(na.omit(lastlist$clep)), width="140px"),
        selectInput("barl_style", "Select Identity or Similarity", choices = c("Identity","Similarity")),
        br(),
        plotOutput("barlev", width = "1500px", height = "1000px"),
        downloadButton("downloadBarLevel")
      )
    })
    
    output$g6 <- renderUI({
      conditionalPanel(
        condition = "input.Barplot % 2 == 1",
        br(),
        br(),
        numericInput("bar_cluster", "Select cluster to show:", 5,  min = 1, max = max(df$clusters), width="140px"),
        selectInput("barcl_style", "Select Identity or Similarity", c("Identity","Similarity")),
        br(),
        plotOutput("barcl",width = "600px", height = "350px"),
        downloadButton("downloadBarCluster")
      )
    })
    
    output$g7 <- renderUI({
      conditionalPanel(
        condition = "input.Amino % 2 == 1",
        br(),
        br(),
        numericInput("amino_level", "Select level to show:", 5,  min = 1, max = max(na.omit(lastlist$clep)), width="140px"),
        br(),
        tableOutput('aminolev')
      )
    })
    
    output$g8 <- renderUI({
      conditionalPanel(
        condition = "input.Amino % 2 == 1",
        br(),
        br(),
        numericInput("amino_cluster", "Select cluster to show:", 5,  min = 0, max = max(df$clusters), width="140px"),
        br(),
        tableOutput('aminocl')
      )
    })
    
    output$g9 <- renderUI({
      conditionalPanel(
        condition = "input.AminoIde % 2 == 1",
        br(),
        tableOutput('pin'),
        br(),
        downloadButton("downloadallEmPin", "Download")
      )
    })
    
    output$g10 <- renderUI({
      conditionalPanel(
        condition = "input.AminoIde % 2 == 1",
        br(),
        br(),
        numericInput("aminoide_level", "Select level to show:", 5,  min = 1, max = max(na.omit(lastlist$clep)), width="140px"),
        selectInput("aminoide_lever", "Select Identity or Similarity", choices = c("Identity","Similarity")),
        br(),
        tableOutput('idenaminolev')
      )
    })
    
    output$g11 <- renderUI({
      conditionalPanel(
        condition = "input.AminoIde % 2 == 1",
        br(),
        br(),
        numericInput("aminoide_cluster", "Select cluster to show:", 5,  min = 1, max = max(df$clusters), width="140px"),
        selectInput("aminoide_clver", "Select Identity or Similarity", choices = c("Identity","Similarity")),
        br(),
        tableOutput('idenaminocl')
      )
    })
    
    
    output$g12 <- renderUI({
      conditionalPanel(
        condition = "input.ConsensusSeq % 2 == 1",
        br(),
        br(),
        numericInput("opt_level", "Select level to show:", 5,  min = 1, max = max(na.omit(lastlist$clep)), width="140px"),
        br("Identity Consensus Sequence for a specific level"),
        textOutput('opt1lev'),
        br("Identity Consensus Sequence with higher frequency letters for a specific level"),
        textOutput('opt1Newlev'),
        br("Similarity Consensus Sequence for a specific level"),
        textOutput('opt2lev'),
        br("Similarity Consensus Sequence with higher frequency groups for a specific level"),
        textOutput('opt2Newlev')
      )
    })
    
    output$g13 <- renderUI({
      conditionalPanel(
        condition = "input.ConsensusSeq % 2 == 1",
        br(),
        br(),
        numericInput("opt_cluster", "Select cluster to show:", 5,  min = 0, max = max(df$clusters), width="140px"),
        br("Identity Consensus Sequence for a specific cluster"),
        textOutput('opt1cl'),
        br("Identity Consensus Sequence with higher frequencye letters for a specific cluster"),
        textOutput('opt1Newcl'),
        br("Similarity Consensus Sequence for a specific cluster"),
        textOutput('opt2cl'),
        br("Similarity Consensus Sequence with higher frequency groups for a specific cluster"),
        textOutput('opt2Newcl')
      )
    })
    
    output$g14 <- renderUI({
      conditionalPanel(
        condition = "input.Identity % 2 == 1",
        selectInput("ident", "Select Identity or Similarity", choices = c("Identity","Similarity")),
        br(),
        plotOutput('idplot'),
        downloadButton("downloadIdentityPlot")
      )
    })
    
    output$g15 <- renderUI({
      conditionalPanel(
        condition = "input.Network % 2 == 1",
        br(),
        br(),
        inline( numericInput("net_level", "Select level to show:", 2,  min = 0, max = max(na.omit(lastlist$clep)), width="140px")),
        inline(selectInput("net_style", "Select Whole or Partial", choices = c("whole","partial"))),
        inline(selectInput("thr_style", "Select Distance, String Identity Similarity, String Group Similarity or Final distance", choices = c("Distance","StrIdentSimilarity","StrGroupSimilarity","FinalDis"))),
        inline(numericInput("thr_value", "Select threshold value:", 2,  min = 1, max = 100, width="140px")),
        inline(textOutput("dis_max")),
        checkboxInput("table_show", "Show Arithmetic Table", FALSE),
        checkboxInput("sil", "Show Network using Silhouette", FALSE),
        #checkboxInput("thr_en", "Show Network using thresholds", FALSE),
        #numericInput("thr_1", "Select cluster to show:", 5,  min = 0, max = 100, width="140px"),
        #numericInput("thr_2", "Select cluster to show:", 10,  min = 0, max = 100, width="140px"),
        #observeEvent(input$thr_en, {
        #  shinyjs::disable("thr_1")
        #  shinyjs::disable("thr_2")
        #}),
        br(),
        actionButton("Apply", "Apply Changes", style="color: #fff; background-color: #179B12; border-color: #fff"),
        plotOutput('net', width = "1500px", height = "1000px"),
        dataTableOutput('net_table')
      )
    })
    
    #session$onSessionEnded(function() {
    #  log.txt <- tic.log(format = TRUE)
    #  write.table(unlist(log.txt), sprintf('name_%s,seqnumber_%d,timestamp_%s.txt', fname, length(udata$AA.JUNCTION), format(Sys.time(),"%Y-%m-%d_%H.%M.%S")), sep="\t")
    #})
    
    })    
 })
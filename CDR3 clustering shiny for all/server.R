


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
      hide("begin")
      udata=read.csv(paste0("data/",input$Dataset), header = TRUE, sep = ";")
      udata$AA.JUNCTION = as.character(udata$AA.JUNCTION)
      udata$clusters = 0 # Initialiaze the column clusters with 0
      udata$level.0 = 0 # Initialize the column of cl.0 with 0
      udata$temp= 0 # Creating a temp column with 0
      mm = 0
      listq = list("ggdf" = ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy,"udata" = udata,"permat"= NA, "persim" = NA, "br" = br, "cl" = cl, "met" = met, "ep" = ep , "clep" = clep, "nn" = nn, "sumper" = NA,"sumper2" = NA, "ela" = NA, "cel" = NA,"endper" = endper, "last" = last)
      lastlist = Matrices(listq,FALSE,let,sim,d,input$algorithm_type,input$column_num)
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
    
      mm = 1
    
      output$message <- renderText({
        if(is.null(mm)) return("sdasad")
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
        BarCl(input$bar_cluster,perlist,persimlist,Clus,let,sim,input$barcl_style,lastlist)
      })
    
      output$downloadBarCluster <- downloadHandler(
        filename = function(){paste0("barPlot_cluster_",input$bar_cluster,".png")},
        content = function(file) {
            png(file,width=1000, height=550)
            BarCl(input$bar_cluster,perlist,persimlist,Clus,let,sim,input$barcl_style,lastlist)
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
        Opt2(AminoLev(input$opt2_level,lastlist,df,Clus),sim,FALSE,altsim)
      })
    
      output$opt2cl <- renderText({
        if (is.null(dir()) | is.null(input$Dataset)) return()
        Opt2(AminoCl(input$opt2_cluster,lastlist,df),sim,FALSE,altsim)
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
        Netw(input$net_level,input$thr_value,input$thr_style,input$net_style,df,lastlist,Clus,sim,altsim)
      })
      
    output$g1 <- renderUI({
      conditionalPanel(
        condition = "input.Tree % 2 == 1",
        br(),
        br(),
        numericInput("tree_level", "Select max level of the tree:", 5,  min = 1, max = max(lastlist$clep), width="140px"),
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
        numericInput("logo_level", "Select level to show:", 5,  min = 1, max = max(lastlist$clep), width="140px"),
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
        numericInput("bar_level", "Select level to show:", 5,  min = 1, max = max(lastlist$clep), width="140px"),
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
        numericInput("amino_level", "Select level to show:", 5,  min = 1, max = max(lastlist$clep), width="140px"),
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
        numericInput("aminoide_level", "Select level to show:", 5,  min = 1, max = max(lastlist$clep), width="140px"),
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
        condition = "input.ComLetters % 2 == 1",
        br(),
        br(),
        numericInput("opt1_level", "Select level to show:", 5,  min = 1, max = max(lastlist$clep), width="140px"),
        br(),
        textOutput('opt1lev')
      )
    })
    
    output$g13 <- renderUI({
      conditionalPanel(
        condition = "input.ComLetters % 2 == 1",
        br(),
        br(),
        numericInput("opt1_cluster", "Select cluster to show:", 5,  min = 0, max = max(df$clusters), width="140px"),
        br(),
        textOutput('opt1cl')
      )
    })
    
    output$g14 <- renderUI({
      conditionalPanel(
        condition = "input.ComGroups % 2 == 1",
        br(),
        br(),
        numericInput("opt2_level", "Select level to show:", 5,  min = 1, max = max(lastlist$clep), width="140px"),
        br(),
        textOutput('opt2lev')
      )
    })
    
    output$g15 <- renderUI({
      conditionalPanel(
        condition = "input.ComGroups % 2 == 1",
        br(),
        br(),
        numericInput("opt2_cluster", "Select cluster to show:", 5,  min = 1, max = max(df$clusters), width="140px"),
        br(),
        textOutput('opt2cl')
      )
    })
    
    output$g16 <- renderUI({
      conditionalPanel(
        condition = "input.Identity % 2 == 1",
        selectInput("ident", "Select Identity or Similarity", choices = c("Identity","Similarity")),
        br(),
        plotOutput('idplot'),
        downloadButton("downloadIdentityPlot")
      )
    })
    
    output$g17 <- renderUI({
      conditionalPanel(
        condition = "input.Network % 2 == 1",
        br(),
        br(),
        numericInput("net_level", "Select level to show:", 2,  min = 0, max = max(lastlist$clep), width="140px"),
        selectInput("net_style", "Select Whole or Partial", choices = c("whole","partial")),
        selectInput("thr_style", "Select Distance, String Similarity or Final distance", choices = c("Distance","StrSimilarity","FinalDis")),
        numericInput("thr_value", "Select threshold value:", 2,  min = 1, max = max(lastlist$clep), width="140px"),
        br(),
        plotOutput('net', width = "1500px", height = "1000px")
      )
    })
    
    })    
 })
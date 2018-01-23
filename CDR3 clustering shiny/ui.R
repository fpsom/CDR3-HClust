#install.packages("shiny")
library(shiny)
library(shinyFiles)
library(shinyjs)
library("shinyBS")
library("DT")

#install.packages("stringr")
#install.packages("dplyr")
#install.packages("entropy")
#install.packages("ggplot2")
#install.packages("ggseqlogo")
#install.packages("gridExtra")
#install.packages("cluster")
#install.packages("seqinr")
#install.packages("collapsibleTree")
#install.packages("data.tree")
#install.packages("DiagrammeR")
#install.packages('stringdist')
#install.packages("igraph")
#install.packages("plsgenomics")
#install.packages("networkD3")
#install.packages('ndtv', dependencies=T)
#install.packages("htmlwidget")
library('networkD3')
library('stringdist')
library("plsgenomics")
library("igraph")
library('network')
library('ndtv')
library("stringr")
library("dplyr")
library("entropy")
library("ggplot2")
library("ggseqlogo")
library("gridExtra")
library("cluster")
library("seqinr")
library("DiagrammeR")
library("collapsibleTree")
library("data.tree")
#library("htmlwidget")
#jsResetCode <- "shinyjs.reset = function() {history.go(0)}"


function(request) {navbarPage(
  #shinyjs::useShinyjs(),
  "CDR3 Clustering",
  id = "navbar",
  position = "fixed-top",  
  inverse=TRUE,
  header = tagList(
    useShinyjs()
    #extendShinyjs("www/app-shinyjs.js", functions = c("updateHistory"))
  ),
  
  tabPanel("Home", value = "home",
           
    tags$head(tags$style(
     HTML('
          #sidebar {
          background-color: #ffffff;
          }

          body, label, input, button, select { 
          font-family: "Arial";
          }')
          )),
  
  tags$head(
    tags$style(HTML("
                    .multicol {
                    -height: 150px;
                    -webkit-column-count: 2; /* Chrome, Safari, Opera */ 
                    -moz-column-count: 2;    /* Firefox */ 
                    column-count: 2; 
                    -moz-column-fill: auto;
                    -column-fill: auto;
                    }
                    "))
    ),
  
  mainPanel( 
    width=9,
    h3("Import Data"),
    br()
  ),
  
  mainPanel(
    width = 9,   
    
    #textOutput("sourced"),
    
    br(),
    
    h4("Select the directory where the folders of the patients' data are located"),
    shinyDirButton("dir", "Choose directory", "Upload"),
    
    #uiOutput("uiInputFiles"),
    
    br(),
    br(),
    uiOutput("uiDatasets"),
    
    br(),
    br(),
    
    actionButton("begin", "Run", 
                 style="color: #fff; background-color: #179B12; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.begin % 2 == 1"
    ),
    textOutput("message"),
    
    br(),
    br(),
    
    actionButton("Tree", "Show Tree", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.Tree % 2 == 1",
      br(),
      br(),
      numericInput("tree_level", "Select max level of the tree:", 5,  min = 1, max = 19, width="140px"),
      br(),
      grVizOutput("tre", width = "1500px", height = "1000px")
    ),

    br(),
    br(),
    
    actionButton("Coltree", "Show Collapsible Tree", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.Coltree % 2 == 1",
      br(),
      br(),
      collapsibleTreeOutput("coltre", width = "1500px", height = "1000px")
    ),
    
    br(),
    br(),
    
    actionButton("Logo", "Show Logo", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.Logo % 2 == 1",
      br(),
      br(),
      numericInput("logo_level", "Select level to show:", 5,  min = 1, max = 19, width="140px"),
      br(),
      plotOutput("logolev"),
      br(),
      downloadButton("downloadLogoLevel")  
    ),
    
    conditionalPanel(
      condition = "input.Logo % 2 == 1",
      br(),
      br(),
      numericInput("logo_cluster", "Select cluster to show:", 5,  min = 1, max = 166, width="140px"),
      br(),
      plotOutput("logocl"),
      br(),
      downloadButton("downloadLogoCluster")
    ),
    
    
    br(),
    br(),
    
    actionButton("Barplot", "Show Barplot", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.Barplot % 2 == 1",
      br(),
      br(),
      numericInput("bar_level", "Select level to show:", 5,  min = 1, max = 18, width="140px"),
      selectInput("barl_style", "Select Identity or Similarity", choices = c("Identity","Similarity")),
      br(),
      plotOutput("barlev", width = "1500px", height = "1000px"),
      downloadButton("downloadBarLevel")
    ),
    
    conditionalPanel(
      condition = "input.Barplot % 2 == 1",
      br(),
      br(),
      numericInput("bar_cluster", "Select cluster to show:", 5,  min = 1, max = 163, width="140px"),
      selectInput("barcl_style", "Select Identity or Similarity", c("Identity","Similarity")),
      br(),
      plotOutput("barcl",width = "600px", height = "350px"),
      downloadButton("downloadBarCluster")
    ),
    
    br(),
    br(),
    
    actionButton("Amino", "Show Amino Elements", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.Amino % 2 == 1",
      br(),
      br(),
      numericInput("amino_level", "Select level to show:", 5,  min = 1, max = 19, width="140px"),
      br(),
      tableOutput('aminolev')
    ),
    
    conditionalPanel(
      condition = "input.Amino % 2 == 1",
      br(),
      br(),
      numericInput("amino_cluster", "Select cluster to show:", 5,  min = 0, max = 166, width="140px"),
      br(),
      tableOutput('aminocl')
    ),
    
    br(),
    br(),
    
    actionButton("AminoIde", "Show Amino Cluster Identity and Similarity", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.AminoIde % 2 == 1",
      br(),
      tableOutput('pin'),
      br(),
      downloadButton("downloadallEmPin", "Download")
    ),
    
    conditionalPanel(
      condition = "input.AminoIde % 2 == 1",
      br(),
      br(),
      numericInput("aminoide_level", "Select level to show:", 5,  min = 1, max = 19, width="140px"),
      selectInput("aminoide_lever", "Select Identity or Similarity", choices = c("Identity","Similarity")),
      br(),
      tableOutput('idenaminolev')
    ),
    
    conditionalPanel(
      condition = "input.AminoIde % 2 == 1",
      br(),
      br(),
      numericInput("aminoide_cluster", "Select cluster to show:", 5,  min = 1, max = 166, width="140px"),
      selectInput("aminoide_clver", "Select Identity or Similarity", choices = c("Identity","Similarity")),
      br(),
      tableOutput('idenaminocl')
    ),
    
    br(),
    br(),
    
    actionButton("ComLetters", "Show Common Letters", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.ComLetters % 2 == 1",
      br(),
      br(),
      numericInput("opt1_level", "Select level to show:", 5,  min = 1, max = 19, width="140px"),
      br(),
      textOutput('opt1lev')
    ),
    
    conditionalPanel(
      condition = "input.ComLetters % 2 == 1",
      br(),
      br(),
      numericInput("opt1_cluster", "Select cluster to show:", 5,  min = 0, max = 166, width="140px"),
      br(),
      textOutput('opt1cl')
    ),
    
    br(),
    br(),
    
    actionButton("ComGroups", "Show Common Similarity Groups", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.ComGroups % 2 == 1",
      br(),
      br(),
      numericInput("opt2_level", "Select level to show:", 5,  min = 1, max = 19, width="140px"),
      br(),
      textOutput('opt2lev')
    ),
    
    conditionalPanel(
      condition = "input.ComGroups % 2 == 1",
      br(),
      br(),
      numericInput("opt2_cluster", "Select cluster to show:", 5,  min = 1, max = 166, width="140px"),
      br(),
      textOutput('opt2cl')
    ),
    
    br(),
    br(),
    
    actionButton("Identity", "Show Identity and Similarity Plot", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.Identity % 2 == 1",
      selectInput("ident", "Select Identity or Similarity", choices = c("Identity","Similarity")),
      br(),
      plotOutput('idplot'),
      downloadButton("downloadIdentityPlot")
    ),
    
    br(),
    br(),
    
    actionButton("Network", "Show Network", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.Network % 2 == 1",
      br(),
      br(),
      numericInput("net_level", "Select level to show:", 2,  min = 0, max = 18, width="140px"),
      selectInput("net_style", "Select Whole or Partial", choices = c("whole","partial")),
      selectInput("thr_style", "Select Distance, String Similarity or Final distance", choices = c("Distance","StrSimilarity","FinalDis")),
      numericInput("thr_value", "Select threshold value:", 2,  min = 1, max = 18, width="140px"),
      br(),
      plotOutput('net', width = "1500px", height = "1000px")
      #downloadButton("downloadBarLevel")
    ),
    
    br(),
    br(),
    
    actionButton("NetworkInt", "Show Interactive Network", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.NetworkInt % 2 == 1",
      br(),
      br(),
      numericInput("netInt_level", "Select level to show:", 2,  min = 0, max = 18, width="140px"),
      selectInput("netInt_style", "Select Whole or Partial", choices = c("whole","partial")),
      selectInput("thrInt_style", "Select Distance, String Similarity or Final distance", choices = c("Distance","StrSimilarity","FinalDis")),
      numericInput("thrInt_value", "Select threshold value:", 2,  min = 1, max = 18, width="140px"),
      br(),
      h2('ndtv-d3 interactive network animation'),
      tags$style(HTML(".tooltip {opacity: 1}")), # stop boostrap css from messing up the tooltip in the widget
      ndtvAnimationWidgetOutput('netInt', width = "1500px", height = "1000px")
    )

  )
  )
  
)}
    
    
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
    
    actionButton("Algorithmtype", "Select Algorithm Type and columns you want to exclude from Algorithm", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    conditionalPanel(
      condition = "input.Algorithmtype % 2 == 1",
      br(),
      br(),
      selectInput("algorithm_type", "Select Identity or Similarity", choices = c("Identity","Similarity")),
      numericInput("column_num", "Select columns to exclude from the algorithm [0 - value]:", 0,  min = 0, max = 20, width="140px"),
      br(),
      
      actionButton("begin", "Run", 
                   style="color: #fff; background-color: #179B12; border-color: #fff"),
      
      conditionalPanel(
        condition = "input.begin % 2 == 1"
      )
    ),
    
    textOutput("message"),
    
    br(),
    
    actionButton("Tree", "Show Tree", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    uiOutput("g1"),
    
    br(),
    
    actionButton("Coltree", "Show Collapsible Tree", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    uiOutput("g2"),
    
    br(),
    
    actionButton("Logo", "Show Logo", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    uiOutput("g3"),
    uiOutput("g4"),
    
    br(),
    
    actionButton("Barplot", "Show Barplot", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    uiOutput("g5"),
    uiOutput("g6"),
    
    br(),
    
    actionButton("Amino", "Show Amino Elements", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    uiOutput("g7"),
    uiOutput("g8"),
    
    br(),
    
    actionButton("AminoIde", "Show Amino Cluster Identity and Similarity", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    
    uiOutput("g9"),
    uiOutput("g10"),
    uiOutput("g11"),
  
    br(),
    
    actionButton("ComLetters", "Show Common Letters", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    uiOutput("g12"),
    uiOutput("g13"),
    
    br(),
    
    actionButton("ComGroups", "Show Common Similarity Groups", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    uiOutput("g14"),
    uiOutput("g15"),
 
    br(),
    
    actionButton("Identity", "Show Identity and Similarity Plot", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    uiOutput("g16"),
    
    br(),
    
    actionButton("Network", "Show Network", 
                 style="color: #fff; background-color: #5F021F; border-color: #fff"),
    uiOutput("g17")

  )
  )
  
)}
    
    
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("stringr","dplyr","entropy","ggplot2","ggseqlogo","gridExtra","cluster","seqinr","collapsibleTree","data.tree","DiagrammeR","stringdist","igraph","networkD3","plsgenomics","shinycssloaders","shiny","shinyFiles","shinyjs","shinyBS","DT","plotly","xtable","tictoc")
ipak(packages)
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
                           style="color: #fff; background-color: #179B12; border-color: #fff")
              
            ),
            
            conditionalPanel(
              condition = "input.begin% 2 == 1",
              withSpinner(textOutput("message"),proxy.height = "50px")),
            
            br()
            
          )
            ),
  
  tabPanel("Tree", value = "TreeTab",
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Tree", "Show Tree", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g1")
           )),
  tabPanel("Coltree", value = "ColtreeTab", 
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Coltree", "Show Collapsible Tree", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g2")
           )),
  tabPanel("Logo", value = "LogoTab",
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Logo", "Show Logo", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g3"),
             uiOutput("g4")
           )),
  tabPanel("Barplot", value = "BarplotTab",
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Barplot", "Show Barplot", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g5"),
             uiOutput("g6")
           )),
  tabPanel( "Amino", value = "AminoTab",
            mainPanel( 
              br(),
              br(),
              br(),
              br(),
              actionButton("Amino", "Show Amino Elements", style="color: #fff; background-color: #5F021F; border-color: #fff"),
              uiOutput("g7"),
              uiOutput("g8")
            )),
  tabPanel("AminoIde", value = "AminoIdeTab", 
           mainPanel(
             br(),
             br(),
             br(),
             br(),
             actionButton("AminoIde", "Show Amino Cluster Identity and Similarity", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g9"),
             uiOutput("g10"),
             uiOutput("g11")
           )),
  tabPanel("ComLetters", value = "ComLettersTab", 
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("ComLetters", "Show Common Letters", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g12"),
             uiOutput("g13")
           )),
  tabPanel("ComGroups", value = "ComGroupsTab", 
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("ComGroups", "Show Common Similarity Groups", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g14"),
             uiOutput("g15")
           )),
  tabPanel("Identity", value = "IdentityTab", 
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Identity", "Show Identity and Similarity Plot", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g16")
           )),
  tabPanel("Network", value = "NetworkTab",
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Network", "Show Network", style="color: #fff; background-color: #5F021F; border-color: #fff"),
            uiOutput("g17")
           ))
  
)}


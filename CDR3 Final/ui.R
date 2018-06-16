######################################## Algorithm Implementation ########################################
# A Multi-metric Algorithm for Hierarchical Clustering of Same-length Protein Sequences - Shiny ui.R
# Tsarouchis Sotrios - Filippos
# email: sotirisftsar@gmail.com
# AEM: 7999
######################################## Install or load necessary libraries ########################################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("stringr","dplyr","entropy","ggplot2","ggseqlogo","gridExtra","cluster","seqinr","collapsibleTree","data.tree","DiagrammeR","stringdist","igraph","networkD3","plsgenomics","shinycssloaders","shiny","shinyFiles","shinyjs","shinyBS","DT","plotly","xtable","tictoc","data.table")
ipak(packages)

######################################## Creation of Tabs ########################################
function(request) {navbarPage(
  #shinyjs::useShinyjs(),
  "CDR3 Clustering",
  id = "navbar",
  position = "fixed-top",  
  inverse=TRUE,
  header = tagList(
    useShinyjs()
  ),
######################################## Home Tab ########################################
  tabPanel("Home", value = "home",
          mainPanel(
            width = 9,   
            fluidPage(theme = "bootstrap.css",
              includeScript("./www/text.js"),
              titlePanel("Folder content upload"),
                fluidRow(
                  column(4,
                    wellPanel(
                      tags$div(class="form-group shiny-input-container", 
                        tags$div(tags$label("File input")),
                        tags$div(tags$label("Choose folder", class="btn btn-primary",
                           tags$input(id = "fileIn", webkitdirectory = TRUE, type = "file", style="display: none;", onchange="pressed()"))),
                        tags$label("No folder choosen", id = "noFile"),
                        tags$div(id="fileIn_progress", class="progress progress-striped active shiny-file-input-progress",
                          tags$div(class="progress-bar")
                        )     
                      ),
                      verbatimTextOutput("results")
                    )
                  ),
                  column(8,
                    checkboxInput("fdet", "Show File Details", FALSE),
                    tabsetPanel(
                      tabPanel("Files table", dataTableOutput("tbl")),
                      tabPanel("Files list", dataTableOutput("tbl2"))
                    )
                    
                    
                  )
                )
            ),
            HTML("<script type='text/javascript' src='getFolders.js'></script>"),
            
            
            br(),
            br(),
            
            actionButton("Algorithmtype", "Select Algorithm Type and columns you want to exclude from Algorithm", 
                         style="color: #fff; background-color: #5F021F; border-color: #fff"),
            
            conditionalPanel(
              condition = "input.Algorithmtype % 2 == 1",
              br(),
              br(),
              tags$div(style="display:inline-block;", selectInput("algorithm_type", "Select Identity or Similarity", choices = c("Identity","Similarity")) ),
              tags$div(style="display:inline-block;",numericInput("column_num", "Select columns to exclude from start [0 - value]:", 0,  min = 0, max = 20, width="140px")),
              tags$div(style="display:inline-block;",numericInput("back_num", "Select columns to exclude from end [0 - value]:", 3,  min = 0, max = 20, width="140px")),
              tags$div(style="display:inline-block;",numericInput("backj6_num", "Select columns to exclude from end if J6 [0 - value]:", 3,  min = 0, max = 20, width="140px")),
              br(),
              tags$div(style="display:inline-block;",checkboxInput("thr_en", "Use thresholds", FALSE)),
              tags$div(style="display:inline-block;",numericInput("thr_1", "Select threshold percentage for number of sequences", 7,  min = 0, max = 100, width="140px")),
              tags$div(style="display:inline-block;",numericInput("thr_2", "Select threshold percentage for Identity or Similarity", 5,  min = 0, max = 100, width="140px")),
              tags$div(style="display:inline-block;",numericInput("thr_level", "Select from which level the threshold will be applied", 0,  min = 0, max = 100, width="140px")),
              
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
  
######################################## Tree Tab ########################################
  tabPanel("Tree", value = "TreeTab",
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Tree", "Show Tree", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g1")
           )
  ),
######################################## Collapsible Tree Tab ########################################
  tabPanel("Coltree", value = "ColtreeTab", 
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Coltree", "Show Collapsible Tree", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g2")
           )
  ),
######################################## Logo Tab ########################################
  tabPanel("Logo", value = "LogoTab",
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Logo", "Show Logo", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g3"),
             uiOutput("g4")
           )
  ),
######################################## Barplot Tab ########################################
  tabPanel("Barplot", value = "BarplotTab",
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Barplot", "Show Barplot", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g5"),
             uiOutput("g6")
           )
  ),
######################################## Amino Tab ########################################
  tabPanel( "Amino", value = "AminoTab",
            mainPanel( 
              br(),
              br(),
              br(),
              br(),
              actionButton("Amino", "Show Amino Elements", style="color: #fff; background-color: #5F021F; border-color: #fff"),
              #uiOutput("g7"),
              uiOutput("g8")
            )
  ),
######################################## Identity Tab ########################################
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
           )
  ),
######################################## Consensus Tab ########################################

  tabPanel("ConsensusSeq", value = "ConsensusSeqTab", 
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("ConsensusSeq", "Show Consensus Sequence", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g12"),
             uiOutput("g13")
           )
  ),
######################################## Identity Plot Tab ########################################
  tabPanel("Identity", value = "IdentityTab", 
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Identity", "Show Identity and Similarity Plot", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             uiOutput("g14")
           )
  ),
######################################## Network Tab ########################################
  tabPanel("Network", value = "NetworkTab",
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("Network", "Show Network", style="color: #fff; background-color: #5F021F; border-color: #fff"),
            uiOutput("g15")
           )
  )
  
)}
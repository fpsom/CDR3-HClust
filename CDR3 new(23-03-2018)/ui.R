ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("stringr","dplyr","entropy","ggplot2","ggseqlogo","gridExtra","cluster","seqinr","collapsibleTree","data.tree","DiagrammeR","stringdist","igraph","networkD3","plsgenomics","shinycssloaders","shiny","shinyFiles","shinyjs","shinyBS","DT","plotly","xtable","tictoc")
ipak(packages)


function(request) {navbarPage(
  #shinyjs::useShinyjs(),
  "CDR3 Clustering",
  id = "navbar",
  position = "fixed-top",  
  inverse=TRUE,
  header = tagList(
    useShinyjs()
  ),
  
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
           )
  ),
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
  tabPanel( "Amino", value = "AminoTab",
            mainPanel( 
              br(),
              br(),
              br(),
              br(),
              actionButton("Amino", "Show Amino Elements", style="color: #fff; background-color: #5F021F; border-color: #fff"),
              uiOutput("g7"),
              uiOutput("g8")
            )
  ),
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
  tabPanel("ConsensusSeq", value = "ConsensusSeqTab", 
           mainPanel( 
             br(),
             br(),
             br(),
             br(),
             actionButton("ConsensusSeq", "Show Consensus Sequence", style="color: #fff; background-color: #5F021F; border-color: #fff"),
             splitLayout( uiOutput("g12"),
                          uiOutput("g13"))
           )
  ),
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
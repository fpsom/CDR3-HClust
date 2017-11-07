# Insert Data
library("stringr")
library("dplyr")
#data <- read.csv(file.choose(), header = TRUE, sep = ";")
data <- read.csv("data/SampleData.csv", header = TRUE, sep = ";")


# Usefull Data
udata <- data[c(1,3)]
udata$AA.JUNCTION <- as.character(udata$AA.JUNCTION)

# Function Lett finds the used letters in our sequences 
Lett <- function(df){
  j = 1
  z = vector(length = 26)
  for (i in 1:26){
    temp = str_count(df$AA.JUNCTION, LETTERS[i])
    if (sum(temp)==0){
      z[j] = i
      j = j+1
    }
  }
  let = LETTERS[-z]
  Matrices(df,let)   # Call the function Matrices
}

# Function Matrices compute the apsolute matrix and the matrix with percentages
Matrices <- function(df1,let){
  mymat <- matrix(0,nrow=length(let), ncol=str_length(df1$AA.JUNCTION[1]))
  permat <- matrix(0,nrow=length(let), ncol=str_length(df1$AA.JUNCTION[1]))
  for (i in 1:length(let)) {
    for (j in 1:str_length(df1$AA.JUNCTION[1])) {
      for (k in 1:length(df1$AA.JUNCTION)) {
        xar <- str_sub(df1$AA.JUNCTION[k],j,j)
        if ( xar == let[i]){
          mymat[i,j] = mymat[i,j] + 1
        }
      }
      permat[i,j] = (mymat[i,j] / length(df1$AA.JUNCTION)) * 100
    }
  }
  Choice(permat,df1,let) # Call the function Choice
}

# Function Choise choose which matrix cell will be used for the division of the data 
Choice <- function(pin,df2,let){
  cel = which(pin == max(pin), arr.ind = TRUE)
  if (max(pin) == 100){ # We exclude the 100 % from the max values
    cel = which(pin == max(pin[pin!=max(pin)]), arr.ind = TRUE) # the desired cell
  }
  Divide(df2,cel,let) # Call the function Divide
}

# Function Divide divide the data into 2 new data frames
Divide <- function(df3,cel,let){
  df4 = filter(df3,str_detect(str_sub(df3$AA.JUNCTION,cel[1,2],cel[1,2]), let[cel[1,1]])) # we use the first max value
  df5 = filter(df3,!str_detect(str_sub(df3$AA.JUNCTION,cel[1,2],cel[1,2]), let[cel[1,1]]))
  Control(df4,df5,cel,let)   # Call the function Control
} 

# Function Control check if the data frame's length is small enough
Control <- function(df6,df7,cel,let){
  if (length(df6$AA.JUNCTION) <2){
    return (list (df6,df7))
  } else if (length(df7$AA.JUNCTION) <2) {
    return (list (df6,df7))
  }
  Matrices(df6,let) # Call the function Matrices
  Matrices(df7,let) # Call the function Matrices
}


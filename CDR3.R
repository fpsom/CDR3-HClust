# Insert Data
library("stringr")
library("dplyr")
library("entropy")
data <- read.csv(file.choose(), header = TRUE, sep = ";")
#data <- read.csv("data/SampleData.csv", header = TRUE, sep = ";")


# Usefull Data
udata <- data[c(1,3)]
udata$AA.JUNCTION <- as.character(udata$AA.JUNCTION)
let = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y") # The letters matrix
br = 0 # initial value of branch
cl = 0 # initial value of new clusters 
udata$clusters = 0 # initialiaze the column clusters with 0
udata$End = FALSE # initialize the column End with FALSE

Matrices <- function(df1,br1,cl1){
  mymat <- matrix(0,nrow=length(let) + 1, ncol=str_length(df1$AA.JUNCTION[1]))
  permat <- matrix(0,nrow=length(let) + 1, ncol=str_length(df1$AA.JUNCTION[1]))
  for (i in 1:str_length(df1[df1$clusters == br1,]$AA.JUNCTION[1])) {
    f <- table(str_sub(df1[df1$clusters == br1,]$AA.JUNCTION,i,i))       # A table with the letters from a specific position 
    k = 1
    for(j in 1:length(let)){
      if (names(f[k]) == let[j] && k <= length(f)){      # We create the mymat with the values of the previous table and zeroes to other letters 
        mymat[j,i] = f[k]
        k = k + 1
      } else{
        mymat[j,i] = 0
      }
      permat[j,i] = (mymat[j,i] / length(df1[df1$clusters == br1,]$AA.JUNCTION)) * 100        # We calculate the percentage matrix
    }
    permat[j+1,i] = entropy(f,base=exp(1))              # We calculate the entropy for every position (row 21)
  }
  Choice(df1,permat,br1,cl1)
}


Choice <- function(df2,pin,br2,cl2){
  cel = which(pin == max(pin), arr.ind = TRUE)
  ela = 1
  if (max(pin) == 100){     # We exclude the 100 % from the max values
    cel = which(pin == max(pin[pin!=max(pin)]), arr.ind = TRUE)    # the desired cell
  } 
  if ((length(cel)/2) > 1){
    for(i in 2:(length(cel)/2)){
      if (pin[21,cel[i,2]] < pin[21,cel[ela,2]]) {          # if the percentage is the same (cel multidimensional) we keep the cel with the lowest entropy
        ela = i
      }
    }
  }
  Divide(df2,br2,cl2,cel,ela)
}


Divide <- function(df3,br3,cl3,cel,ela){      # We change the clusters column for the 2 new Clusters
  df3[df3$clusters == br3,]$clusters <- ifelse(str_detect(str_sub(df3[df3$clusters == br3,]$AA.JUNCTION,cel[ela,2],cel[ela,2]), let[cel[ela,1]]), cl3+1 ,cl3+2)
  cl3 = cl3 + 2                               # the value of the new clusters
  Control(df3,br3,cl3)
}


Control <- function(df4,br4,cl4){
  if (nrow(df4[df4$clusters == cl4,]) <= 2){      # We check if the number of sequences for this Cluster is <= 2
    df4[df4$clusters == cl4,]$End = TRUE          # we change the End column
    br4 = cl4 - 1                                 # the new branch
    if (nrow(df4[df4$clusters == cl4-1,]) <= 2){  # We check if the number of sequences for the previous Cluster is <= 2
      df4[df4$clusters == cl4-1,]$End = TRUE      # we change the End column
      if (nrow(df4[df4$End == FALSE,]) <= 2 ){    # We check if the number of sequences with the End column FALSE (remaining sequences) is <= 2
        return(df4)                               # We end our algortihm
      }
      br4 = max((df4[df4$End == FALSE,])$clusters)   # the new branch
    }
  } else {
    br4 = cl4   # the new branch
  }
  Matrices(df4,br4,cl4)
}


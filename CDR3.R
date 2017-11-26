# Insert Data
library("stringr")
library("dplyr")
library("entropy")
#data <- read.csv(file.choose(), header = TRUE, sep = ";")
data <- read.csv("data/SampleData.csv", header = TRUE, sep = ";")


# Usefull Data
udata <- data[c(1,3)]
udata$AA.JUNCTION <- as.character(udata$AA.JUNCTION)
# A list with the similarity groups
sim = list("F","W",c("A","I","L","V"),c("M","C"),"P","G","Y",c("T","S"),c("H","K","R"),c("E","D"),c("Q","N"))
# Naming the group of similarities
names(sim) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
# A table with the letters
let = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y") # The letters matrix
# An empty vector wich will store the level for every cluster
clep = vector('numeric')
br = 0  # initial value of branch
cl = 0  # initial value of new clusters
met = 0 # initial value of level capacity counter
ep = 1  # initial value of level
udata$clusters = 0 # initialiaze the column clusters with 0
udata$cl.0 = 0 # initialize the column of cl.0 with 0
udata$temp= 0 # creating a temp column with 0
d = 0

# Function Matrices compute the apsolute matrix, the matrix with percentages, the absolute similarity matrix and the percentage similarity matrix
Matrices <- function(){
  mymat <<- matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
  permat <<- matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
  simmat <<- matrix(0,nrow = length(sim),ncol =str_length(udata$AA.JUNCTION[1]))
  persim <<- matrix(0,nrow = length(sim)+1,ncol =str_length(udata$AA.JUNCTION[1]))
  for (i in 1:str_length(udata[udata$clusters == br,]$AA.JUNCTION[1])) {
    f <- table(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,i,i)) # A table with the letters from a specific position
    k = 1
    sum =0
    for(j in 1:length(let)){
      if (names(f[k]) == let[j] && k <= length(f)){ # We create the mymat with the values of the previous table and zeroes to other letters
        mymat[j,i] <<- f[k]
        d <<- str_which(sim,let[j]) # We determine the group of every letter and we create the simmat
        simmat[d,i] <<- simmat[d,i] + f[k]
        k = k + 1
      } else{
        mymat[j,i] <<- 0
      }
      permat[j,i] <<- (mymat[j,i] / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100   # We calculate the percentage matrix
      persim[d,i] <<- (simmat[d,i] / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100  # We calculate the percentage matrix for similarity groups
    }
    permat[j+1,i] <<- entropy(f,base=exp(1))                      # We calculate the entropy for every position (row 21)
    persim[length(sim)+1,i] <<- entropy(simmat[,i],base = exp(1)) # We calculate the entropy for every position (row 12)
  }
  Choice()
}

# Function Choice choose which matrix cell will be used for the division of the data 
Choice <- function(){
  cel <<- which(permat == max(permat), arr.ind = TRUE)
  ela <<- 1
  poss = max(permat) 
  if (max(permat) == 100){ # We exclude the 100 % from the max values
    cel <<- which(permat == max(permat[permat!=max(permat)]), arr.ind = TRUE) # the desired cell
    poss = max(permat[permat!=max(permat)])
  }
  if ((length(cel)/2) > 1){
    eqcol = vector('numeric') # A vector for the cell's rows with the same entropy
    j = 1
    for(i in 1:(length(cel)/2)){
      if (permat[21,cel[i,2]] < permat[21,cel[ela,2]]) {  # if the percentage is the same (cel multidimensional) we keep the cel with the lowest entropy
        ela <<- i
        eqcol[] = 0
        eqcol[1] = i
        j = 2
      }else if(permat[21,cel[i,2]] == permat[21,cel[ela,2]]){
        eqcol[j] = i
        j = j+1
      }
    }
    eqcol = eqcol[-which(eqcol == 0)]
    if (length(eqcol)>1){ # if the vector has 2 or more numbers means that we have columns with the same entropy
      ela <<- eqcol[1]
      meg = persim[str_which(sim,let[cel[ela,1]]),cel[ela,2]]
      for(i in 2:length(eqcol)){
        if(persim[str_which(sim,let[cel[eqcol[i],1]]),cel[eqcol[i],2]] > meg){ # We find the maximum percentage of the similarity percentage marix
          meg = persim[str_which(sim,let[cel[eqcol[i],1]]),cel[eqcol[i],2]]
          ela <<- eqcol[i]
        } else if (persim[str_which(sim,let[cel[eqcol[i],1]]),cel[eqcol[i],2]] == meg){ # If the percentage is the same we use the entropy
          if(persim[12,cel[eqcol[i],2]] < persim[12,cel[ela,2]]){
            ela <<- i
          }
        }
      }
      
    }
  }
  Divide()
}

# Function Divide divide the data into 2 new clusters creating 2 columns with name (cl.cluster_number) and updating the column Clusters
Divide <- function(){
  udata$temp <<- udata[,str_which(names(udata),sprintf('cl.%d', br))]
  # Creating the new column with the appropriate data
  udata[udata$temp == br,]$temp <<- ifelse(str_detect(str_sub(udata[udata$temp == br,]$AA.JUNCTION,cel[ela,2],cel[ela,2]), let[cel[ela,1]]), cl+1 ,br)
  # Naming the column as cl.cluster_number
  names(udata)[length(udata)] <<- sprintf('cl.%d', cl+1)
  udata$temp <<- udata[,str_which(names(udata),sprintf('cl.%d', br))]
  udata[udata$temp == br,]$temp <<- ifelse(str_detect(str_sub(udata[udata$temp == br,]$AA.JUNCTION,cel[ela,2],cel[ela,2]), let[cel[ela,1]]), br ,cl+2)
  names(udata)[length(udata)] <<- sprintf('cl.%d', cl+2)
  # Updating the cluster column
  udata[udata$clusters == br,]$clusters <<- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,cel[ela,2],cel[ela,2]), let[cel[ela,1]]), cl+1 ,cl+2)
  Control()
}

# Function Control check if the data frame's length is small enough
Control <- function(){
  clep[cl+1] <<- ep # Level of the cluster cl+1  
  clep[cl+2] <<- ep # Level of the cluster cl+2
  br <<- br + 1 # Increase the branch by 1
  cl <<- cl + 2 # Increase the cluster by 2
  met <<- met + 2 # Increase the counter by 2
  while (length(which(udata$clusters == br)) <= 2 && br != cl){ # While the number of sequences in the branch is less than 2, go to the next branch and change counter 
    if( ((clep[br] < clep[br+1]) && (sum(clep == clep[br]) != (2^clep[br]))) == TRUE ){ # If the next branch is in the next level
      met <<- 0
      ep <<- ep + 1
      br <<- br + 1
    }else{
      br <<- br + 1
      met <<- met +2 
    }
  }
  if( br+1> cl && br !=0){ # When the branch reach the last value of cluster return the data
    return(udata)
  }
  if( met == geomSeq(1,2,1,50)[ep+1]){ # When the counter reaches the end value (geometric sequence) we increase the level counter
    met <<- 0
    ep <<- ep + 1
  }
  Matrices()
}

# A function that generates a geometric sequence
geomSeq <- function(start,ratio,begin,end){
  begin=begin-1
  end=end-1
  start*ratio**(begin:end)
}

# A function which visualize the sequences of a data frame using common letters (i.e. "A _ _ _ _ K R _ _ _ Q _ Y Y Y _ _ _ T _")
Opt <- function(df){
  xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
  for (i in 1:str_length(df$AA.JUNCTION[1])) {
    f <- table(str_sub(df$AA.JUNCTION,i,i))
    xar[i] = "_"
    if (f[1] == nrow(df)){
      xar[i] = names(f[1])
    }
    
  }
  paste(xar,collapse = ' ')
}

# A function which visualize the sequences of a data frame using similarity groups (i.e. "A _ _ _ _ K Am _ Ba _ Ac _ Y Y Y _ _ _ T _")
Opt2 <- function(df){
  xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
  for (i in 1:str_length(df$AA.JUNCTION[1])) {
    f <- table(str_sub(df$AA.JUNCTION,i,i))
    xar[i] = "_"
    if (f[1] == nrow(df)){
      xar[i] = names(f[1])
    }else{
      y = TRUE
      d = str_which(sim,names(f[1]))
      for(j in 2:length(f)){
        y = y && str_detect(sim[d],names(f[j]))
      }
      if( y == TRUE){
        xar[i] = names(sim[d])
      } 
    }
    
  }
  paste(xar,collapse = ' ')
}

# Using the Opt and Opt2 for the data in cluster 33 for example
x= as.data.frame(udata[udata$cl.33 == 33,]$AA.JUNCTION)
names(x)[1]= "AA.JUNCTION"

# Creating a plot for cluster 34 for example
ggseqlogo(udata[udata$cl.34 == 34,]$AA.JUNCTION)

# An alternative plot for a cluster (We must determine the cluster's permat)
barplot(permat, col=heat.colors(length(permat[,1])-1), width=2)
legend("bottomright",inset=c(0,0), fill=heat.colors(length(permat[,1])-1), legend=let)

# Insert Data

#install.packages("stringr")
#install.packages("dplyr")
#install.packages("entropy")
#install.packages("ggplot2")
#install.packages("ggseqlogo")
#install.packages("gridExtra")
#install.packages("cluster")
#install.packages("seqinr")
library("stringr")
library("dplyr")
library("entropy")
library("ggplot2")
library("ggseqlogo")
library("gridExtra")
library("cluster")
library("seqinr")
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
br = 0  # Initial value of branch
cl = 0  # Initial value of new clusters
met = 0 # Initial value of level capacity counter
ep = 1  # Initial value of level
udata$clusters = 0 # Initialiaze the column clusters with 0
udata$level.0 = 0 # Initialize the column of cl.0 with 0
udata$temp= 0 # Creating a temp column with 0
d = 0
nn = FALSE # Initial value for the condition sumper < endper%
endper = 90 # Set the percentage, which ends the programm 
listxx = list() # Initialize a list for saving the permat of all branches
dfsum = data.frame(sumper = numeric(0),sumper2 = numeric(0),branch = numeric(0), len = numeric(0)) 
# Initalize the first list
listq = list("dfsum" = dfsum,"list" = listxx,"udata" = udata,"permat"= NA, "persim" = NA, "br" = br, "cl" = cl, "met" = met, "ep" = ep , "clep" = clep, "nn" = nn, "sumper" = NA,"sumper2" = NA, "ela" = NA, "cel" = NA,"endper" = endper)


# Function Matrices compute the apsolute matrix, the matrix with percentages, the absolute similarity matrix and the percentage similarity matrix
Matrices <- function(list1,leaf){
  udata = list1$udata
  permat = list1$permat
  persim = list1$persim
  br = list1$br
  cl = list1$cl
  listxx = list1$list
  mymat = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
  permat = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
  simmat = matrix(0,nrow = length(sim),ncol =str_length(udata$AA.JUNCTION[1]))
  persim = matrix(0,nrow = length(sim)+1,ncol =str_length(udata$AA.JUNCTION[1]))
  rownames(mymat) = c(let,"Entropy")
  rownames(permat) = c(let,"Entropy")
  rownames(simmat) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
  rownames(persim) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am","Entropy")
  for (i in 1:str_length(udata[udata$clusters == br,]$AA.JUNCTION[1])) {
    f <- table(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,i,i)) # A table with the letters from a specific position
    k = 1
    sum =0
    for(j in 1:length(let)){
      if (names(f[k]) == let[j] && k <= length(f)){ # We create the mymat with the values of the previous table and zeroes to other letters
        mymat[j,i] = f[k]
        d = str_which(sim,let[j]) # We determine the group of every letter and we create the simmat
        simmat[d,i] = simmat[d,i] + f[k]
        k = k + 1
      } else{
        mymat[j,i] = 0
      }
      permat[j,i] = (mymat[j,i] / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100   # We calculate the percentage matrix
      persim[d,i] = (simmat[d,i] / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100  # We calculate the percentage matrix for similarity groups
    }
    permat[j+1,i] = entropy(f,base=exp(1))                      # We calculate the entropy for every position (row 21)
    persim[length(sim)+1,i] = entropy(simmat[,i],base = exp(1)) # We calculate the entropy for every position (row 12)
  }
  listxx$temp = permat
  names(listxx)[length(listxx)] = sprintf('permat_br.%d', br) # Save the permat with this format
  result1=list("dfsum" = list1$dfsum,"list" = listxx, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep" = list1$ep, "clep" = list1$clep, "nn" = list1$nn, "sumper" =list1$sumper,"sumper2" = list1$sumper2, "ela" = list1$ela, "cel" = list1$cel, "endper" = list1$endper)
  if ( leaf == TRUE){
    Finish(result1,TRUE)
  }else{
    Finish(result1,FALSE)
  }
}


Finish <- function(list2,leaf){
  udata = list2$udata
  permat = list2$permat
  persim = list2$persim
  sumper = list2$sumper # The total percentage of permat 
  sumper2 = list2$sumper2 # The total percentage of persim
  endper = list2$endper
  dfsum = list2$dfsum
  br = list2$br
  t1 =which(permat[-length(permat),] == 100,arr.ind = TRUE)
  sumper = (length(as.numeric(t1[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
  #sumper2 = NA
  t2 =which(persim[-length(persim),] == 100,arr.ind = TRUE)
  sumper2 = (length(as.numeric(t2[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
  
  if (sumper > endper){
    nn = TRUE  # When nn = TRUE the percentage of sumper < endper%
    #t2 =which(persim[-length(persim),] == 100,arr.ind = TRUE)
    #sumper2 = (length(as.numeric(t2[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
    if (sumper2 > endper){
      return(list2)  # End of the programm, returning the final list
    }
  }
  vv = length(udata[udata$clusters == br,]$AA.JUNCTION)
  dfsum[nrow(dfsum) + 1,] = c(sumper,sumper2,br,vv)
  result2 = list("dfsum" = dfsum,"list" = list2$list,"udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = list2$br, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper)
  if ( leaf == TRUE){
    nn = FALSE
    result2 = list("dfsum" = dfsum,"list" = list2$list,"udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = list2$br, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper)
    return(result2)
  }else{
    Choice(result2)
  }
}


# Function Choice choose which matrix cell will be used for the division of the data 
Choice <- function(list3){
  nn = list3$nn
  if(nn == TRUE){
    permat = list3$persim # If sumper < endper% we want to check only the persim matrix
  }else{
    permat =list3$permat  # Else permat and if it is necessary the persim matrix
    persim = list3$persim
  }
  cel = which(permat == max(permat), arr.ind = TRUE)
  ela = 1
  poss = max(permat) 
  if (max(permat) == 100){ # We exclude the 100 % from the max values
    cel = which(permat == max(permat[permat!=max(permat)]), arr.ind = TRUE) # The desired cell
    poss = max(permat[permat!=max(permat)])
  }
  if ((length(cel)/2) > 1){
    eqcol = vector('numeric') # A vector for the cell's rows with the same entropy
    j = 1
    for(i in 1:(length(cel)/2)){
      if (permat[21,cel[i,2]] < permat[21,cel[ela,2]]) {  # If the percentage is the same (cel multidimensional) we keep the cel with the lowest entropy
        ela = i
        eqcol[] = 0
        eqcol[1] = i
        j = 2
      }else if(permat[21,cel[i,2]] == permat[21,cel[ela,2]]){
        eqcol[j] = i
        j = j+1
      }
    }
    eqcol <- eqcol[-which(eqcol == 0)]
    if (length(eqcol)>1 && nn == FALSE){ # If the vector has 2 or more numbers means that we have columns with the same entropy and nn = FALSE in order not to double check the persim
      ela = eqcol[1]
      meg = persim[str_which(sim,let[cel[ela,1]]),cel[ela,2]]
      for(i in 2:length(eqcol)){
        if(persim[str_which(sim,let[cel[eqcol[i],1]]),cel[eqcol[i],2]] > meg){ # We find the maximum percentage of the similarity percentage marix
          meg = persim[str_which(sim,let[cel[eqcol[i],1]]),cel[eqcol[i],2]]
          ela = eqcol[i]
        } else if (persim[str_which(sim,let[cel[eqcol[i],1]]),cel[eqcol[i],2]] == meg){ # If the percentage is the same we use the entropy
          if(persim[12,cel[eqcol[i],2]] < persim[12,cel[ela,2]]){
            ela = i
          }
        }
      }
      
    }
  }
  nn = FALSE # Return nn in it's original value 
  result3 = list("dfsum" = list3$dfsum,"list" = list3$list,"udata" = list3$udata,"permat"= list3$permat, "persim" = list3$persim, "br" = list3$br, "cl" = list3$cl, "met" = list3$met, "ep"= list3$ep, "clep" = list3$clep,"nn" = nn, "sumper" = list3$sumper, "sumper2" = list3$sumper2, "ela" = ela, "cel" = cel,"endper" = list3$endper)
  Divide(result3)
}


# Function Divide divide the data into 2 new clusters and updates the column with the right level 
Divide <- function(list4){
  udata = list4$udata
  br = list4$br
  cel = list4$cel
  ela = list4$ela
  cl = list4$cl
  met = list4$met
  ep = list4$ep
  if (met == 0){ # If we need a new level, then we create a new column with its name (level.ep)
    udata$temp = NA
    names(udata)[length(udata)] = sprintf('level.%d', ep)
  }
  # The update of the right column
  x1 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,cel[ela,2],cel[ela,2]), let[cel[ela,1]]), "TRUE")
  y1 = udata[udata$clusters == br,]$AA.JUNCTION
  z1 = y1[x1]
  x2 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,cel[ela,2],cel[ela,2]), let[cel[ela,1]]), "FALSE")
  y2 = udata[udata$clusters == br,]$AA.JUNCTION
  z2 = y2[x2]
  for(i in 1:length(z1)){
    udata[str_which(udata$AA.JUNCTION,z1[i]),sprintf('level.%d', ep)] = cl+1
  }
  for(i in 1:length(z2)){
    udata[str_which(udata$AA.JUNCTION,z2[i]),sprintf('level.%d', ep)] = cl+2
  }
  # Updating the cluster column
  udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,cel[ela,2],cel[ela,2]), let[cel[ela,1]]), cl+1 ,cl+2) 
  result4 = list("dfsum" = list4$dfsum,"list" = list4$list,"udata" = udata,"permat"= list4$permat, "persim" = list4$persim, "br" = list4$br, "cl" = list4$cl, "met" = list4$met, "ep"= list4$ep, "clep" = list4$clep,"nn" = list4$nn, "sumper" = list4$sumper, "sumper2" = list4$sumper2, "ela" = list4$ela, "cel" = list4$cel,"endper" = list4$endper)
  Control(result4)
}


# Function Control update the branch, the counter and the cluster numbers in order to begin a new valid sub-division
Control <- function(list5){
  udata = list5$udata
  br = list5$br
  cel = list5$cel
  ela = list5$ela
  cl = list5$cl
  met = list5$met
  clep = list5$clep
  ep = list5$ep
  clep[cl+1] = ep # Level of the cluster cl+1  
  clep[cl+2] = ep # Level of the cluster cl+2
  br = br + 1 # Increase the branch by 1
  cl = cl + 2 # Increase the cluster by 2
  met = met + 2 # Increase the counter by 2
  list5 = list("dfsum" = list5$dfsum,"list" = list5$list,"udata" = list5$udata,"permat"= list5$permat, "persim" = list5$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list5$nn, "sumper" = list5$sumper, "sumper2" = list5$sumper2, "ela" = list5$ela, "cel" = list5$cel,"endper" = list5$endper)
  while (length(which(udata$clusters == br)) <= 2){ # While the number of sequences in the branch is less than 2, go to the next branch and change counter 
    list5 = Matrices(list5, TRUE)
    if( ((clep[br] < clep[br+1]) && (sum(clep == clep[br]) != (2^clep[br]))) == TRUE ){ # If the next branch is in the next level
      met = 0
      ep = ep + 1
      br = br + 1
    }else{
      br = br + 1
      met = met +2 
    }
    list5 = list("dfsum" = list5$dfsum,"list" = list5$list,"udata" = list5$udata,"permat"= list5$permat, "persim" = list5$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list5$nn, "sumper" = list5$sumper, "sumper2" = list5$sumper2, "ela" = list5$ela, "cel" = list5$cel,"endper" = list5$endper)
  }
  if( met == geomSeq(1,2,1,50)[ep+1]){ # When the counter reaches the end value (geometric sequence) we increase the level counter
    met = 0
    ep = ep + 1
  }
  result5 = list("dfsum" = list5$dfsum,"list" = list5$list,"udata" = list5$udata,"permat"= list5$permat, "persim" = list5$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list5$nn, "sumper" = list5$sumper, "sumper2" = list5$sumper2, "ela" = list5$ela, "cel" = list5$cel,"endper" = list5$endper)
  Matrices(result5,FALSE)
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

# A command for beggining our programm 
lastlist = Matrices(listq,FALSE)
# The final name of udata data frame 
df = lastlist$udata
# A list with the permat matrix for every branch
perlist = lastlist$list

# Create custom colour scheme
cs1 = make_col_scheme(chars=c("F","W","A","I","L","V","M","C","P","G","Y","T","S","H","K","R","E","D","Q","N"),
                      cols=c("#1E90FF", "#BA55D3", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#C6E2FF", "#C6E2FF", "#FFD700", "#00EE00", "#C1FFC1", "#54FF9F", "#54FF9F", "#FF0000", "#FF0000", "#FF0000", "#FFD700", "#FFD700", "#ED9121", "#ED9121"))
# Creating a plot for cluster 5 for example
ggseqlogo(df[df$level.2 == 5,]$AA.JUNCTION, method = "prob", col_scheme=cs1)

# A function to plot with logo level lev
Logo <- function(lev){
  t1 = which(lastlist$clep == lev)
  listff = list()
  if(length(t1) %% 3 == 0){
    nc = length(t1)%/%3
  }else{
    nc = length(t1)%/%3 + 1
  }
  
  for(i in 1:length(t1)){
    x1 = as.data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", lev))] == t1[i], ]$AA.JUNCTION))
    names(x1)[1]= "AA.JUNCTION"
    x1 = as.character(x1$AA.JUNCTION)
    listff$temp = x1
    names(listff)[length(listff)] = sprintf('Cluster.%d', t1[i])
  }
  ggseqlogo(listff, ncol=nc, method = "prob",col_scheme=cs1)
}

# Create the grid (we need to miximize the pop-up window in order to have a perfect fit)
get( getOption( "device" ) )()
# Creating a plot for level 5 for example
Logo(6)

# An alternative plot for cluster 3 (We must determine the cluster's permat)
par(mar=c(3,3,4,4),xpd=TRUE)
output <- matrix(unlist(perlist[3]), ncol = 20, byrow = FALSE)
barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d',3)) 
legend("topright",inset=c(-0.03,-0.23), fill=heat.colors(length(output[,1])-1), legend=let,cex = 0.6)


# A function to plot with barplot level lev
Bar <- function(lev){
  t2 = which(lastlist$clep == lev)
  if(length(t2) %% 3 == 0){
    par(mfrow = c(length(t2)%/%3,3))
  }else{
    par(mfrow = c(length(t2)%/%3 + 1,3))
  }
  for(i in 1:length(t2)){
    ar = str_which(names(perlist),as.character(t2[i]))[1]
    par(xpd=TRUE)
    output <- matrix(unlist(perlist[ar]), ncol = 20, byrow = FALSE)
    barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d', t2[i])) 
    legend("topright",inset=c(-0.03,-0.23), fill=heat.colors(length(output[,1])-1), legend=let,cex = 0.6)
  }
}

# Create the grid (we need to miximize the pop-up window in order to have a perfect fit)
get( getOption( "device" ) )()
# Plot level 5
Bar(5)


# Sequences and Id's for specific cluster
clust = 15
lastlist$clep[clust]
x3 = data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$Sequence.ID), na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$AA.JUNCTION))
names(x3) = c("Sequence.ID","AA.JUNCTION")
se = x3$Sequence.ID
aa = as.character(x3$AA.JUNCTION)
aa = as.list(aa)
names(aa) = se
x3 # Print in console
write.fasta(aa, names= names(aa) ,file = sprintf("Cluster.%d", clust)) # Create a fasta file
Opt(x3) 
opt2(x3)

# Sequences and Id's for specific level
level = 7
x4 = na.omit(as.matrix(df[which((names(df) == sprintf("level.%d", level)))]))
x4 = df[x4,][1:2]
se2 = x4$Sequence.ID
aa2 = as.character(x4$AA.JUNCTION)
aa2 = as.list(aa2)
names(aa2) = se2
x4 # Print in console
write.fasta(aa2, names= names(aa2) ,file = sprintf("Level.%d", level)) # Create a fasta file
Opt(x4)
opt2(x4)

# Plot the identity(sumper)
ff = lastlist$dfsum
ff$level[1]= 0
ff$level[2:(length(lastlist$clep[ff$branch]) +1 )] = lastlist$clep[ff$branch]
par(xpd=TRUE)
plot(ff$branch,ff$sumper,col = ff$level, pch = ff$level, main = "Identity Plot", xlab = "Number of Clusters", ylab = "Identity")
legend("topright",inset=c(-0.05,0.11), col = unique(ff$level), pch = unique(ff$level), legend= unique(sprintf("Level.%d", ff$level)),cex = 0.6)
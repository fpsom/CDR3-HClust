#install.packages(c("stringr","dplyr","entropy","ggplot2","ggseqlogo","gridExtra","cluster","seqinr","collapsibleTree","data.tree","DiagrammeR"))
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
#install.packages('stringdist')
#install.packages("igraph")
#install.packages("plsgenomics")
library('stringdist')
library("plsgenomics")
library("igraph")
library("stringr")
library("dplyr")
library("entropy")
library("ggplot2")
library("ggseqlogo")
library("gridExtra")
library("cluster")
library("seqinr")
library("collapsibleTree")
library("data.tree")
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
listyy = list()
dfsum = data.frame(sumper = numeric(0),sumper2 = numeric(0),branch = numeric(0), len = numeric(0)) 
# Initalize the first list
listq = list("dfsum" = dfsum,"list" = listxx, "listn" = listyy,"udata" = udata,"permat"= NA, "persim" = NA, "br" = br, "cl" = cl, "met" = met, "ep" = ep , "clep" = clep, "nn" = nn, "sumper" = NA,"sumper2" = NA, "ela" = NA, "cel" = NA,"endper" = endper)


# Function Matrices compute the apsolute matrix, the matrix with percentages, the absolute similarity matrix and the percentage similarity matrix
Matrices <- function(list1,leaf){
  udata = list1$udata
  permat = list1$permat
  persim = list1$persim
  br = list1$br
  cl = list1$cl
  listxx = list1$list
  listyy = list1$listn
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
  listyy$temp = persim
  names(listyy)[length(listyy)] = sprintf('persim_br.%d', br)
  result1=list("dfsum" = list1$dfsum,"list" = listxx,"listn" = listyy, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep" = list1$ep, "clep" = list1$clep, "nn" = list1$nn, "sumper" =list1$sumper,"sumper2" = list1$sumper2, "ela" = list1$ela, "cel" = list1$cel, "endper" = list1$endper)
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
  result2 = list("dfsum" = dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = list2$br, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper)
  if ( leaf == TRUE){
    nn = FALSE
    result2 = list("dfsum" = dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = list2$br, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper)
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
      if (permat[nrow(permat),cel[i,2]] < permat[nrow(permat),cel[ela,2]]) {  # If the percentage is the same (cel multidimensional) we keep the cel with the lowest entropy
        ela = i
        eqcol[] = 0
        eqcol[1] = i
        j = 2
      }else if(permat[nrow(permat),cel[i,2]] == permat[nrow(permat),cel[ela,2]]){
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
          if(persim[nrow(persim),cel[eqcol[i],2]] < persim[nrow(persim),cel[ela,2]]){
            ela = i
          }
        }
      }
      
    }
  }
  nn = FALSE # Return nn in it's original value 
  result3 = list("dfsum" = list3$dfsum,"list" = list3$list, "listn" = list3$listn, "udata" = list3$udata,"permat"= list3$permat, "persim" = list3$persim, "br" = list3$br, "cl" = list3$cl, "met" = list3$met, "ep"= list3$ep, "clep" = list3$clep,"nn" = nn, "sumper" = list3$sumper, "sumper2" = list3$sumper2, "ela" = ela, "cel" = cel,"endper" = list3$endper)
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
  result4 = list("dfsum" = list4$dfsum,"list" = list4$list, "listn" = list4$listn, "udata" = udata,"permat"= list4$permat, "persim" = list4$persim, "br" = list4$br, "cl" = list4$cl, "met" = list4$met, "ep"= list4$ep, "clep" = list4$clep,"nn" = list4$nn, "sumper" = list4$sumper, "sumper2" = list4$sumper2, "ela" = list4$ela, "cel" = list4$cel,"endper" = list4$endper)
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
  list5 = list("dfsum" = list5$dfsum,"list" = list5$list, "listn" = list5$listn, "udata" = list5$udata,"permat"= list5$permat, "persim" = list5$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list5$nn, "sumper" = list5$sumper, "sumper2" = list5$sumper2, "ela" = list5$ela, "cel" = list5$cel,"endper" = list5$endper)
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
    list5 = list("dfsum" = list5$dfsum,"list" = list5$list, "listn" = list5$listn, "udata" = list5$udata,"permat"= list5$permat, "persim" = list5$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list5$nn, "sumper" = list5$sumper, "sumper2" = list5$sumper2, "ela" = list5$ela, "cel" = list5$cel,"endper" = list5$endper)
  }
  if( met == geomSeq(1,2,1,50)[ep+1]){ # When the counter reaches the end value (geometric sequence) we increase the level counter
    met = 0
    ep = ep + 1
  }
  result5 = list("dfsum" = list5$dfsum,"list" = list5$list, "listn" = list5$listn, "udata" = list5$udata,"permat"= list5$permat, "persim" = list5$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list5$nn, "sumper" = list5$sumper, "sumper2" = list5$sumper2, "ela" = list5$ela, "cel" = list5$cel,"endper" = list5$endper)
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

persimlist = lastlist$listn
# Create a dataframe with all clusters and their identity and similarity percentage 
ff = lastlist$dfsum
ff$level[1]= 0
ff$level[2:(length(lastlist$clep[ff$branch]) +1 )] = lastlist$clep[ff$branch]  #without leaves


Clus = as.data.frame(matrix(100, ncol = 3, nrow = max(df$clusters)+1))
names(Clus) = c("ClusterId","Identity","Similarity")
Clus$ClusterId = 0:max(df$clusters)
Clus$seqnum = 1
Clus$level = c(0,lastlist$clep)
for(i in 1:length(ff$branch) ){
  ll = which(Clus$ClusterId == ff$branch[i])
  Clus$Identity[ll] = ff$sumper[i]
  Clus$Similarity[ll] = ff$sumper2[i]
  Clus$seqnum[ll] = ff$len[i]
}

AminoCl <- function(clust){
  if(clust == 0){
    x3 = df[1:2]
  }else{
    lastlist$clep[clust]
    x3 = data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$Sequence.ID), na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$AA.JUNCTION))
    names(x3) = c("Sequence.ID","AA.JUNCTION")
    se = x3$Sequence.ID
    aa = as.character(x3$AA.JUNCTION)
    aa = as.list(aa)
    names(aa) = se 
  }
  x3 # Print in console
}

lev = 3 #determine the level

df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
df_args <- c(df[str_which(names(df), "level.")], sep="/")
if(lev == 19){
  df$pathString<- do.call(paste, df_args)
  kk = df$pathString
  for(i in 1:length(kk)){
    temp = str_locate(kk[i],"/NA")
    if(is.na(temp[1]) == FALSE){
      temp2 = str_sub(kk[i], 1, temp[1]-1);
      kk[i] = temp2 
    }
  } 
}else{
  df$pathString<- do.call(paste, df_args)
  #df$pathString <- paste(df$level.0, df$level.1, df$level.2 ,df$level.3 ,df$level.4 ,df$level.5 ,df$level.6 ,df$level.7 ,df$level.8 ,df$level.9,df$level.10 ,df$level.11 ,df$level.12 ,df$level.13 ,df$level.14 ,df$level.15 ,df$level.16, df$level.17, df$level.18, df$level.19 , sep = "/")
  kk = df$pathString
  gg =as.data.frame(str_locate_all(kk,"/"))
  tem = 1
  for (i in 1:length(kk)){
    kk[i] = str_sub(kk[i],1,gg[,tem][lev+1]-1)
    tem = tem + 2
  }
  for(i in 1:length(kk)){
    temp = str_locate(kk[i],"/NA")
    if(is.na(temp[1]) == FALSE){
      temp2 = str_sub(kk[i], 1, temp[1]-1);
      kk[i] = temp2 
    }
  } 
}

df$pathString = kk
x <- ToDataFrameTree(df, "pathstring")
xN <- as.Node(x)

bo = which(lastlist$clep == lev)
ffg = matrix(NA,nrow = max(bo)+1,ncol = (max(bo)+1)) # the distance between 2 Nodes
ffg2 = matrix(NA,nrow = max(bo)+1,ncol = (max(bo)+1)) # the string distance between 2 Nodes (string from function Opt())
ffg3 = matrix(NA,nrow = max(bo)+1,ncol = (max(bo)+1)) # the final distance combining ffg and ffg2

# anw trigwnikos
for(i in 1:(max(bo)+1)){
  tempN1 = FindNode(xN,(sprintf("%d", i-1)))
  tempN1 = tempN1$path
  temp2N1 = Opt(AminoCl(i-1))
  for(j in i:(max(bo)+1)){          # max(df$clusters)+1
    tempN2 = FindNode(xN,(sprintf("%d", j-1)))
    tempN2 = tempN2$path
    tem = which(tempN1 == tempN2)
    tempN1[tem[length(tem)]]
    d1 = length(tempN1) - tem[length(tem)]
    d2 = length(tempN2) - tem[length(tem)]
    ffg[i,j]= d1 + d2
    temp2N2 = Opt(AminoCl(j-1))
    temp2N2 = str_replace_all(temp2N2,"_"," ")
    ffg2[i,j] = stringdist(temp2N1,temp2N2,method = "lv" )
    ffg3[i,j] = -0.25 * ffg[i,j] * 0.7692308 + 0.25 * ffg2[i,j]
  }
}

# Damerau–Levenshtein
matrix.heatmap(ffg) # heatmap for distance between 2 Nodes
matrix.heatmap(ffg2) # heatmap for string distance between 2 Nodes (string from function Opt())
matrix.heatmap(ffg3) # heatmap for final distance combining ffg and ffg2
jjj = matrix(1,nrow = max(bo)+1,ncol = (max(bo)+1))
net0 <- graph_from_adjacency_matrix(jjj,mode = "upper")

bb = vector(length = (max(bo)+1))
bb[1]=0
bb[2:length(bb)] = lastlist$clep[1:max(bo)]
#XRWMA
# Generate colors based on media type:
colrs <- c("#1E90FF", "#BA55D3", "#0000FF", "#557fd2", "#54d17e", "#8aad62", "#C6E2FF", "#e5e234", "#FFD700", "#00EE00", "#C1FFC1", "#ea8509", "#54FF9F", "#FF0000", "#ed3b1c", "#ed1c7a", "#0c0c0c", "#b8d8af", "#ED9121","#45f713")
colrs = colrs[1:(max(bb)+1)]
V(net0)$color <- colrs[bb+1]

# Compute node degrees (#links) and use that to set node size:
deg <- Clus$seqnum[1:length(bb)] / 5
V(net0)$size <- deg

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net0)$label <- NA

# Set edge width based on weight:
hhh = na.omit(as.vector(t(ffg3)))
hhh[which(na.omit(as.vector(t(ffg3))) == 0)] = 0.01
E(net0)$width <- hhh #na.omit(as.vector(t(ffg3)))

#change arrow size and edge color:
E(net0)$arrow.size <- .2
E(net0)$edge.color <- "gray80"

edge.start <- ends(net0, es=E(net0), names=F)[,1]
edge.col <- V(net0)$color[edge.start]
plot(net0, edge.color=edge.col, edge.curved=.1) #plot the network graph
legend(x=-2.5, y=1, colnames(df[str_which(names(df), "level.")])[1:(max(bb)+1)], pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

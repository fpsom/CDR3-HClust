ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("stringr","dplyr","entropy","ggplot2","ggseqlogo","gridExtra","cluster","seqinr","collapsibleTree","data.tree","DiagrammeR","stringdist","igraph","networkD3","plsgenomics","shinycssloaders","shiny","shinyFiles","shinyjs","shinyBS","DT","plotly","xtable","tictoc","data.table","pryr")
ipak(packages)

data <- read.csv(file.choose(), header = TRUE, sep = ";")
udata <- data[c(1,3)]
udata = data[,c("Sequence.ID","J.GENE.and.allele","AA.JUNCTION","V.GENE.and.allele")]
udata$AA.JUNCTION <- as.character(udata$AA.JUNCTION)
sim = list("F","W",c("A","I","L","V"),c("M","C"),"P","G","Y",c("T","S"),c("H","K","R"),c("E","D"),c("Q","N"))
# Naming the group of similarities
names(sim) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
altsim = list("F","W",c("A","I","L","V"),c("M","C"),"P","G","Y",c("T","S"),c("H","K","R"),c("E","D"),c("Q","N"))
# Naming the group of similarities
names(altsim) = c("f","w","a","s","p","g","y","h","b","c","m")
# A table with the letters
let = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y") # The letters matrix
# An empty vector wich will store the level for every cluster
clep = vector('numeric')
br = 0  # Initial value of branch
cl = 0  # Initial value of new clusters
met = 0 # Initial value of level capacity counter
ep = 1  # Initial value of level
d = 0
nn = FALSE # Initial value for the condition sumper < endper%
endper = 90 # Set the percentage, which ends the programm 
listxx = list() # Initialize a list for saving the permat of all branches
listyy = list()
dfsum = data.frame(sumper = numeric(0),sumper2 = numeric(0),branch = numeric(0), len = numeric(0))
dfsd = data.frame(ff = character(0),Average_Identity_Value = numeric(0),Identity_Standar_Deviation = numeric(0),Average_Similarity_Value = numeric(0), Similarity_Standar_Deviation = numeric(0), stringsAsFactors = FALSE) 
ggdf = data.frame(branch = numeric(0), len = numeric(0))
last = 0
flagtic = TRUE
udata$AA.JUNCTION = as.character(udata$AA.JUNCTION)
udata$clusters = 0 # Initialiaze the column clusters with 0
udata$level.0 = 0 # Initialize the column of cl.0 with 0
udata$temp= 0 # Creating a temp column with 0
progend = FALSE
listax = list()
listax$temp = 1:length(udata$AA.JUNCTION)
names(listax)[length(listax)] = sprintf('cl.%d', 0) 
listq = list("listax" = listax,"ggdf" = ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy,"udata" = udata,"permat"= NA, "persim" = NA, "br" = br, "cl" = cl, "met" = met, "ep" = ep , "clep" = clep, "nn" = nn, "sumper" = NA,"sumper2" = NA, "ela" = NA, "cel" = NA,"endper" = endper, "last" = last,"progend" = progend, "leaf" = FALSE,"leaf2" = FALSE)
algo = "Identity"
#algo = "Similarity"
algocol = 2
backcol = 3
backcolj6 = 3
listb = listq
levcut = 3
levcut =levcut +1
leaf2 = FALSE
enthr = FALSE


#list1 = listq
#list1 = result2
# Function Matrices compute the apsolute matrix, the matrix with percentages, the absolute similarity matrix and the percentage similarity matrix
Matrices <- function(list1){
  tic()
  udata = list1$udata
  br = list1$br
  permat = list1$permat
  persim = list1$persim
  cl = list1$cl
  listxx = list1$list
  listyy = list1$listn
  sumper = list1$sumper # The total percentage of permat 
  sumper2 = list1$sumper2 # The total percentage of persim
  endper = list1$endper
  dfsum = list1$dfsum
  nn = list1$nn
  last = list1$last
  progend  = list1$progend
  leaf = list1$leaf
  leaf2 = list1$leaf2

  ind1 = str_which(udata[udata$clusters == br,]$J.GENE.and.allele,"J6")
  qw = 1:length(udata[udata$clusters == br,]$AA.JUNCTION)
  if(length(ind1) == 0){
    ind2 = qw
  }else{
    ind2 = qw[-ind1]
  }

  mymat1 = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
  simmat1 = matrix(0,nrow = length(sim),ncol =str_length(udata$AA.JUNCTION[1]))
  rownames(mymat1) = let
  rownames(simmat1) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
  
  
  permat = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
  rownames(permat) = c(let,"Entropy")
  persim = matrix(0,nrow = length(sim)+1,ncol =str_length(udata$AA.JUNCTION[1]))
  rownames(persim) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am","Entropy")
  print(udata)
  
  if(length(ind1) != 0 ){
    trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION[ind1],"")
    align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
    
    
    for(i in (algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)){
      temptab = plyr::count(align[i],vars = colnames(align)[i])
      names(temptab)[1] = "X1"
      mymat1[which(is.na(match(names(mymat1[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
      #permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
      temptab1 = temptab
      temptab1[1] = names(sim)[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim,x)))]
      temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
      simmat1[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim) == x)),i] = as.vector(unlist(temptab1[2]))
      #persim[length(sim)+1,i] = entropy(simmat1[,i],base = exp(1))
    }
    #permat1[1:length(let),] = (mymat1 / length(udata[udata$clusters == br,]$AA.JUNCTION[ind1])) * 100
    #persim1[1:length(sim),] = (simmat1 / length(udata[udata$clusters == br,]$AA.JUNCTION[ind1])) * 100 
  }
  
  mymat2 = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
  #permat2 = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
  simmat2 = matrix(0,nrow = length(sim),ncol =str_length(udata$AA.JUNCTION[1]))
  #persim2 = matrix(0,nrow = length(sim)+1,ncol =str_length(udata$AA.JUNCTION[1]))
  rownames(mymat2) = let
  #rownames(permat2) = c(let,"Entropy")
  rownames(simmat2) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
  #rownames(persim2) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am","Entropy")
  #print(udata)
  
  if(length(ind2) != 0 ){
    trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION[ind2],"")
    align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
    
    
    for(i in (algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)){
      temptab = plyr::count(align[i],vars = colnames(align)[i])
      names(temptab)[1] = "X1"
      mymat2[which(is.na(match(names(mymat2[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
      permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
      temptab1 = temptab
      temptab1[1] = names(sim)[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim,x)))]
      temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
      simmat2[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim) == x)),i] = as.vector(unlist(temptab1[2]))
      persim[length(sim)+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
    }
    #permat2[1:length(let),] = (mymat2 / length(udata[udata$clusters == br,]$AA.JUNCTION[ind1])) * 100
    #persim2[1:length(sim),] = (simmat2 / length(udata[udata$clusters == br,]$AA.JUNCTION[ind1])) * 100 
  }
  
  mymat = mymat1 + mymat2
  simmat = simmat1 + simmat2
  meg = length(ind1) + length(ind2)
  if(length(ind2) == 0 || backcol == backcolj6){
    permat[1:length(let),(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (mymat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
    persim[1:length(sim),(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (simmat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
  }else{
    permat[1:length(let),(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (mymat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
    persim[1:length(sim),(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (simmat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
    if(backcolj6 != 0){
      permat[1:length(let),(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] = (mymat[,(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] / length(ind2)) * 100
      persim[1:length(sim),(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] = (simmat[,(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] / length(ind2)) * 100
    }
  }
  
  trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION,"")
  align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
  
  for(i in 1:(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]))){
    temptab = plyr::count(align[i],vars = colnames(align)[i])
    names(temptab)[1] = "X1"
    #mymat1[which(is.na(match(names(mymat1[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
    permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
    temptab1 = temptab
    temptab1[1] = names(sim)[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim,x)))]
    temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
    #simmat1[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim) == x)),i] = as.vector(unlist(temptab1[2]))
    persim[length(sim)+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
  }
  
  
  if(algocol != 0 && backcol!=0){
    mymat3 = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
    #permat2 = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
    simmat3 = matrix(0,nrow = length(sim),ncol =str_length(udata$AA.JUNCTION[1]))
    #persim2 = matrix(0,nrow = length(sim)+1,ncol =str_length(udata$AA.JUNCTION[1]))
    rownames(mymat3) = let
    #rownames(permat2) = c(let,"Entropy")
    rownames(simmat3) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
    permat3 = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
    rownames(permat3) = c(let,"Entropy")
    persim3 = matrix(0,nrow = length(sim)+1,ncol =str_length(udata$AA.JUNCTION[1]))
    rownames(persim3) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am","Entropy")
    trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION,"")
    align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
    exc = 1:algocol
    exc2 = (str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol):str_length(udata[udata$clusters == br,]$AA.JUNCTION[1])
    exc = c(exc,exc2)
    for(i in 1:length(exc)){
      temptab = plyr::count(align[exc[i]],vars = colnames(align)[exc[i]])
      names(temptab)[1] = "X1"
      mymat3[which(is.na(match(names(mymat3[,exc[i]]),as.vector(unlist(temptab[1])))) == FALSE),exc[i]] = as.vector(unlist(temptab[2]))
      permat3[length(let) + 1,exc[i]] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
      temptab1 = temptab
      temptab1[1] = names(sim)[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim,x)))]
      temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
      simmat3[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim) == x)),exc[i]] = as.vector(unlist(temptab1[2]))
      persim3[length(sim)+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
    }
    permat3[1:length(let),] = (mymat3 / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100
    persim3[1:length(sim),] = (simmat3 / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100
    permat[1:length(let),exc] = permat3[1:length(let),exc] 
    persim[1:length(sim),exc] = persim3[1:length(sim),exc] 
  }
  
  listxx$temp = permat
  names(listxx)[length(listxx)] = sprintf('permat_br.%d', br) # Save the permat with this format
  listyy$temp = persim
  names(listyy)[length(listyy)] = sprintf('persim_br.%d', br)
  
  ################################ Finish ###########################
  t1 = which(permat[-nrow(permat),] == 100,arr.ind = TRUE)
  sumper = (length(as.numeric(t1[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
  t2 = which(persim[-nrow(persim),] == 100,arr.ind = TRUE)
  sumper2 = (length(as.numeric(t2[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
  vv = length(udata[udata$clusters == br,]$AA.JUNCTION)
  dfsum[nrow(dfsum) + 1,] = c(sumper,sumper2,br,vv)
  if(algo == "Identity"){
    if (sumper > endper && leaf == FALSE){
      nn = TRUE  # When nn = TRUE the percentage of sumper < endper%
      if (sumper2 > endper && leaf == FALSE){
        #progend = TRUE
        print("@@@@@@@@exit 1")
        #list1 = list("listax" = list1$listax,"ggdf" = list1$ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep"= list1$ep, "clep" = list1$clep,"nn" = list1$nn, "sumper" = list1$sumper, "sumper2" = list1$sumper2, "ela" = list1$ela, "cel" = list1$cel,"endper" = list1$endper, "last" = last,"progend" = list1$progend, "leaf" = leaf)

        #en = toc(quiet = TRUE)
        #cat(paste0("Matrices","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
        #return(list1)  # End of the programm, returning the final list
		    leaf2 = TRUE			
      }
    }
  }else{
    if (sumper2 > endper){
      #progend = TRUE
      #list1 = list("listax" = list1$listax,"ggdf" = list1$ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep"= list1$ep, "clep" = list1$clep,"nn" = list1$nn, "sumper" = list1$sumper, "sumper2" = list1$sumper2, "ela" = list1$ela, "cel" = list1$cel,"endper" = list1$endper, "last" = last,"progend" = list1$progend, "leaf" = leaf)

      #en = toc(quiet = TRUE)
      #cat(paste0("Matrices","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
      #return(list1)   # End of the programm, returning the final list
	    leaf2 = TRUE	
	}
  }
  
  result1 = list("listax" = list1$listax,"ggdf" = list1$ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep"= list1$ep, "clep" = list1$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "ela" = list1$ela, "cel" = list1$cel,"endper" = list1$endper, "last" = last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)

  en = toc(quiet = TRUE)
  cat(paste0("Matrices","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
  if ( leaf == TRUE){
    nn = FALSE
    result1 = list("listax" = list1$listax,"ggdf" = list1$ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep"= list1$ep, "clep" = list1$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "ela" = list1$ela, "cel" = list1$cel,"endper" = list1$endper, "last" = last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
    return(result1)
  }else{
    return(result1)
  }
}


#list2 = result1
#list2 = lista
# Function Choice choose which matrix cell will be used for the division of the data 
Choice <- function(list2){
  tic()
  udata = list2$udata
  br = list2$br
  cl = list2$cl
  nn = list2$nn
  cel = list2$cel
  ela = list2$ela
  met = list2$met
  ep = list2$ep
  clep = list2$clep
  ggdf = list2$ggdf
  progend  = list2$progend
  leaf = list2$leaf
  listax = list2$listax
  leaf2 = list2$leaf2
  
  
    if(nn == TRUE || (algo == "Similarity") ){
      permat = list2$persim # If sumper < endper% we want to check only the persim matrix
    }else{
      permat =list2$permat  # Else permat and if it is necessary the persim matrix
      persim = list2$persim
    }
    cel = which(permat[,(algocol + 1):(ncol(permat)-backcol)] == max(permat[,(algocol + 1):(ncol(permat)-backcol)]), arr.ind = TRUE)
    ela = 1
    poss = max(permat[,(algocol + 1):(ncol(permat)-backcol)]) 
    if (max(permat[,(algocol + 1):(ncol(permat)-backcol)]) == 100){ # We exclude the 100 % from the max values
      cel = which(permat[,(algocol + 1):(ncol(permat)-backcol)] == max(permat[,(algocol + 1):(ncol(permat)-backcol)][permat[,(algocol + 1):(ncol(permat)-backcol)]!=max(permat[,(algocol + 1):(ncol(permat)-backcol)])]), arr.ind = TRUE) # The desired cell
      poss = max(permat[,(algocol + 1):(ncol(permat)-backcol)][permat[,(algocol + 1):(ncol(permat)-backcol)]!=max(permat[,(algocol + 1):(ncol(permat)-backcol)])])
    }
  
    #if(poss == 0){
    #  clmax = cl
    #  brtemp = br
    #  progend = TRUE
    #  print(paste0("branch=",brtemp," cl=",cl))
    #  list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
    #  for (i in (br+1):clmax) {
    #    brtemp = i
    #    list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
    #    list2 = Matrices(list2)
    #}
    #print("@@@@@@@@@@@@@@@@@@@@@@@@@@exit 2")
    #en = toc(quiet = TRUE)
    #cat(paste0("Choice","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
    #return(list2)
    #}
if(leaf2 == FALSE && poss != 0){  
  ela = 1
  if(algo == "Identity"){
    if ((length(cel)/2) > 1){
      dddff = min(permat[,(algocol + 1):(ncol(permat)-backcol)][nrow(permat[,(algocol + 1):(ncol(permat)-backcol)]),cel[,2]])
      dddff2 = which(permat[,(algocol + 1):(ncol(permat)-backcol)][nrow(permat[,(algocol + 1):(ncol(permat)-backcol)]),cel[,2]] == dddff)
      ela = dddff2[1]
      if (length(dddff2)>1 && nn == FALSE){ # If the vector has 2 or more numbers means that we have columns with the same entropy and nn = FALSE in order not to double check the persim
        dddff3 = persim[,(algocol + 1):(ncol(permat)-backcol)][sapply(1:length(dddff2), function (x){str_which(sim,let[cel[x,1]])}),cel[,2]]
        dddff3 = max(diag(dddff3))
        dddff4 = which(diag( persim[,(algocol + 1):(ncol(permat)-backcol)][sapply(1:length(dddff2), function (x){str_which(sim,let[cel[x,1]])}),cel[,2]]) == dddff3)
        ela = dddff4[1]
        if(length(dddff4)> 1){
          dddff5 = min(persim[,(algocol + 1):(ncol(permat)-backcol)][nrow(persim[,(algocol + 1):(ncol(permat)-backcol)]),cel[dddff4,2]])
          dddff6 = which(persim[,(algocol + 1):(ncol(permat)-backcol)][nrow(persim[,(algocol + 1):(ncol(permat)-backcol)]),cel[dddff4,2]] == dddff5)
          ela = dddff6[1] 
		   
        }
      }
    }
  }else{
    if ((length(cel)/2) > 1){
      dddff = min(permat[,(algocol + 1):(ncol(permat)-backcol)][nrow(permat[,(algocol + 1):(ncol(permat)-backcol)]),cel[,2]])
      dddff2 = which(permat[,(algocol + 1):(ncol(permat)-backcol)][nrow(permat[,(algocol + 1):(ncol(permat)-backcol)]),cel[,2]] == dddff)
      ela = dddff2[1]
    }
  }
  nn = FALSE # Return nn in it's original value 

  ##################################################### Divide ############################
  if (met == 0){ # If we need a new level, then we create a new column with its name (level.ep)
    udata$temp = NA
    names(udata)[length(udata)] = sprintf('level.%d', ep)
  }

  if (algo == "Identity"){
    x1 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), let[cel[ela,1]]), "TRUE")
    mk1 = length(x1)
    cltp1 = cl + 1
    y1 = udata[udata$clusters == br,]$AA.JUNCTION
    z1 = y1[x1]
    x2 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), let[cel[ela,1]]), "FALSE")
    mk2 = length(x2)
    cltp2 = cl + 2
    y2 = udata[udata$clusters == br,]$AA.JUNCTION
    z2 = y2[x2]
    lengdif = FALSE
    if(mk1 < mk2){
      lengdif = TRUE
      templeng = z1
      z1 = z2
      z2 = templeng
      temp2 = x1
      x1 = x2
      x2 = temp2
      ggdf[nrow(ggdf) + 1,] = c(cltp1,mk2)
      ggdf[nrow(ggdf) + 1,] = c(cltp2,mk1)
    }else{
      ggdf[nrow(ggdf) + 1,] = c(cltp1,mk1)
      ggdf[nrow(ggdf) + 1,] = c(cltp2,mk2)
    }
												  
    listax$temp = as.vector(unlist(listax[sprintf('cl.%d', br)]))[x1]
    names(listax)[length(listax)] = sprintf('cl.%d', cl+1) # Save the permat with this format
    listax$temp = as.vector(unlist(listax[sprintf('cl.%d', br)]))[x2]
    names(listax)[length(listax)] = sprintf('cl.%d', cl+2) # Save the permat with this format
    udata[as.vector(unlist(listax[sprintf('cl.%d', br)]))[x1],sprintf('level.%d', ep)] = cl+1
    udata[as.vector(unlist(listax[sprintf('cl.%d', br)]))[x2],sprintf('level.%d', ep)] = cl+2
    
    if(lengdif == TRUE){
      udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), let[cel[ela,1]]), cl+2 ,cl+1)
    }else{
      udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), let[cel[ela,1]]), cl+1 ,cl+2)
    }
    print(br)
    print(udata$clusters)
  }else{
    strings.to.find = unlist(sim[cel[ela,1]])
    x1 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), str_c(strings.to.find, collapse="|")), "TRUE")
    mk1 = length(x1)
    cltp1 = cl + 1
    y1 = udata[udata$clusters == br,]$AA.JUNCTION
    z1 = y1[x1]
    x2 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), str_c(strings.to.find, collapse="|")), "FALSE")
    mk2 = length(x2)
    cltp2 = cl + 2
    y2 = udata[udata$clusters == br,]$AA.JUNCTION
    z2 = y2[x2]
    lengdif = FALSE
    if(mk1 < mk2){
      lengdif = TRUE
      templeng = z1
      z1 = z2
      z2 = templeng
      temp2 = x1
      x1 = x2
      x2 = temp2
      ggdf[nrow(ggdf) + 1,] = c(cltp1,mk2)
      ggdf[nrow(ggdf) + 1,] = c(cltp2,mk1)
    }else{
      ggdf[nrow(ggdf) + 1,] = c(cltp1,mk1)
      ggdf[nrow(ggdf) + 1,] = c(cltp2,mk2)
    }
    
    listax$temp = as.vector(unlist(listax[sprintf('cl.%d', br)]))[x1]
    names(listax)[length(listax)] = sprintf('cl.%d', cl+1) # Save the permat with this format
    listax$temp = as.vector(unlist(listax[sprintf('cl.%d', br)]))[x2]
    names(listax)[length(listax)] = sprintf('cl.%d', cl+2) # Save the permat with this format
    udata[as.vector(unlist(listax[sprintf('cl.%d', br)]))[x1],sprintf('level.%d', ep)] = cl+1
    udata[as.vector(unlist(listax[sprintf('cl.%d', br)]))[x2],sprintf('level.%d', ep)] = cl+2
    
    if(lengdif == TRUE){
      udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), str_c(strings.to.find, collapse="|")), cl+2 ,cl+1) 
    }else{
      udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), str_c(strings.to.find, collapse="|")), cl+1 ,cl+2) 
    }
  }
    clep[cl+1] = ep # Level of the cluster cl+1  
    clep[cl+2] = ep # Level of the cluster cl+2
    cl = cl + 2 # Increase the cluster by 2
  }else{
    leaf2 = FALSE
    nn = FALSE
  }
  
  ############################################ Control #############################################
  br = br + 1 # Increase the branch by 1
  met = met + 2 # Increase the counter by 2
  list2 = list("listax" = listax,"ggdf" = ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = ela, "cel" = cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
  if( ((clep[br-1] < clep[br]) && (sum(clep == clep[br]) != (2^clep[br]))) == TRUE && (length(which(udata$clusters == br)) > 1) ){ # If the next branch is in the next level
    met = geomSeq(1,2,1,1000)[ep+1]
  }
  while (length(which(udata$clusters == br)) <= 1){ # While the number of sequences in the branch is less than 2, go to the next branch and change counter 
    if(is.na(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]))){
      print("exittttttttt@@@@@@@22")
      #print(paste0("br = ",br))
      clmax = cl
      brtemp = br
      progend = TRUE
      list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
      if(brtemp < clmax){
        for (i in (br+1):clmax) {
          #i = i+1
          brtemp = i
          #print(brtemp)
          #print(clmax)
          list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
          list2 = Matrices(list2)
        }
      }
      en = toc(quiet = TRUE)
      cat(paste0("Choice","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
      return(list2)
    }
    leaf = TRUE
    list2 = list("listax" = listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = ela, "cel" = cel,"endper" = list2$endper, "last" = list2$last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
    list2 = Matrices(list2)
    print(paste0("br",br))
    print(list2$dfsum[br+1,])
    if(is.na(((clep[br] < clep[br+1]) && (sum(clep == clep[br]) != (2^clep[br]))))){
      print("exittttttttt@@@@@@@22")
      #print(paste0("br = ",br))
      clmax = cl
      brtemp = br
      progend = TRUE
      list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
      if(brtemp < clmax){
        for (i in (br+1):clmax) {
          #i = i+1
          brtemp = i
          #print(brtemp)
          #print(clmax)
          list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
          list2 = Matrices(list2)
        }
      }
      en = toc(quiet = TRUE)
      cat(paste0("Choice","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
      return(list2)
    }
    if( ((clep[br] < clep[br+1]) && (sum(clep == clep[br]) != (2^clep[br]))) == TRUE ){ # If the next branch is in the next level
      met = 0
      ep = ep + 1
      br = br + 1
    }else{
      br = br + 1
      met = met +2 
    }
    list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
  }
  if( met == geomSeq(1,2,1,1000)[ep+1]){ # When the counter reaches the end value (geometric sequence) we increase the level counter
    met = 0
    ep = ep + 1
  }
  leaf = FALSE
  result2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = ela, "cel" = cel,"endper" = list2$endper, "last" = list2$last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)

  en = toc(quiet = TRUE)
  cat(paste0("Choice","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
  return(result2)
  
}



# A function that generates a geometric sequence
geomSeq <- function(start,ratio,begin,end){
  begin=begin-1
  end=end-1
  start*ratio**(begin:end)
}

#logfile
logFile = paste0(getwd(),"/log_file ",trunc(as.numeric(Sys.time())),".txt")
cat(paste0("Function","\t","Num of input rows","\t","Num of input columns","\t","Duration","\t","Memory used"), file=logFile, append=FALSE, sep = "\n")

tic()
while (progend == FALSE){
  lista = Matrices(listb)
  progend = lista$progend
  if(progend == TRUE){
    break
  }
  listb = Choice(lista)
  progend = listb$progend
}
en = toc(quiet = TRUE)
cat(paste0("lastlist","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic,"\t",(pryr::mem_used() / 10^6)), file=logFile, append=TRUE, sep = "\n")
lastlist = listb
lastlist2 = listb
lastlist3 = listb
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
Clus$seqnum[1] = nrow(df) 
Clus$seqnum[2:(max(df$clusters)+1)] = lastlist$ggdf$len
Clus$level = c(0,lastlist$clep)
for(i in 1:length(ff$branch) ){
  ll = which(Clus$ClusterId == ff$branch[i])
  Clus$Identity[ll] = ff$sumper[i]
  Clus$Similarity[ll] = ff$sumper2[i]
}


threshold1 = 5
threshold2 = 5
lev = max(na.omit(lastlist$clep))
df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
df_args <- c(df[str_which(names(df), "level.")], sep="/")
if(lev == max(na.omit(lastlist$clep))){
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
if(enthr == TRUE){
  
  matches <- regmatches(unique(x$pathString), gregexpr("[[:digit:]]+", unique(x$pathString)))
  arxpath = unique(x$pathString)
  for(i in 1:length(matches)){
    ggg = as.numeric(unlist(matches[i])) +1
    if(levcut <= (length(ggg)-1)){
      for(j in levcut:(length(ggg)-1)){
        if(algo == "Identity"){
          ddd1 = ((Clus$seqnum[ggg[j]] - Clus$seqnum[ggg[j+1]]) / Clus$seqnum[ggg[j]]) * 100
          ddd2 = (Clus$Identity[ggg[j+1]] - Clus$Identity[ggg[j]]) * 0.5 + (Clus$Similarity[ggg[j+1]] - Clus$Similarity[ggg[j]]) * 0.5
          if(ddd1 < threshold1 || ddd2 < threshold2){
            deik = which(x$pathString == arxpath[i])
            tem = str_locate(arxpath[i],as.character(ggg[j]-1))
            fftemp = str_sub(x$pathString[deik[1]],1,tem[2])
            ff2 = str_sub(x$pathString[deik[1]],tem[1],tem[2])
            ff3 = as.numeric(ff2)
            x$clusters[deik] = ff3
            ff4 = lastlist$clep[ff3]
            ff5 = str_which(names(x),sprintf("level.%d",ff4))
            ff6 = str_which(names(x), "level.")
            if(is.na(ff5[1]) == FALSE){
              x[deik,(ff5[1]+1):max(ff6)] = NA
            }
            x$pathString[deik] = fftemp
            break()
          }
        }else if(algo == "Similarity"){
          ddd1 = ((Clus$seqnum[ggg[j]] - Clus$seqnum[ggg[j+1]]) / Clus$seqnum[ggg[j]]) * 100
          ddd2 = Clus$Similarity[ggg[j+1]] - Clus$Similarity[ggg[j]]
          print(paste0("ddd1=",ddd1," ddd2=",ddd2))
          print(paste0("i=",i," j=",j) )
          if(ddd1 < threshold1 || ddd2 < threshold2){
            deik = which(x$pathString == arxpath[i])
            tem = str_locate(arxpath[i],as.character(ggg[j]-1))
            fftemp = str_sub(x$pathString[deik[1]],1,tem[2])
            ff2 = str_sub(x$pathString[deik[1]],tem[1],tem[2])
            ff3 = as.numeric(ff2)
            x$clusters[deik] = ff3
            ff4 = lastlist$clep[ff3]
            ff5 = str_which(names(x),sprintf("level.%d",ff4))
            ff6 = str_which(names(x), "level.")
            if(is.na(ff5[1]) == FALSE){
              x[deik,(ff5[1]+1):max(ff6)] = NA
            }
            x$pathString[deik] = fftemp
            break()
          }
        }
      }
    }
  }
  xN <- as.Node(x)
  bo = which(lastlist$clep == lev)
  levtel = lastlist$clep
  j = 1
  jfjf = vector()
  for(i in 1:max(bo)){
    temp1 = as.numeric(names(FindNode(xN,(sprintf("%d",i)))$children))
    if(length(temp1) == 1){
      if(temp1[1] %% 2 == 1){
        ff7 = temp1[1] + 1
      }else{
        ff7 = temp1[1] - 1
      }
      ff8 = lastlist$clep[ff7]
      ff9 = str_which(names(x),sprintf("level.%d",ff8))
      deik2 = which(df[ff9] == ff7)
      x[deik2,ff9] = df[deik2,ff9]
      tem = str_locate(df$pathString[deik2[1]],as.character(ff7))
      fftemp = str_sub(df$pathString[deik2[1]],1,tem[2])
      x$pathString[deik2] = fftemp
      jfjf[j] = ff7
      j = j+1
    }
    if(is.null(FindNode(xN,(sprintf("%d",i)))$isLeaf) && is.element(i, jfjf) == FALSE){
      levtel[i] = NA
      Clus[i+1,] = NA
      ff[i+1,] = NA
    }
  }
  xN <- as.Node(x)
}
#df = x
#lastlist$clep = levtel



lev = 4
Den <- function(lev){
  if(flagtic == TRUE) tic(sprintf("Den -- Level: %d ", lev))
  df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
  df_args <- c(df[str_which(names(df), "level.")], sep="/")
  if(lev == max(na.omit(lastlist$clep))){
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
  #par = str_count(df$pathString,"/")
  #str_locate(df$pathString,"/")
  #par = rev(gregexpr("\\/", df$pathString))
  if(enthr == TRUE){
    matches <- regmatches(unique(x$pathString), gregexpr("[[:digit:]]+", unique(x$pathString)))
    arxpath = unique(x$pathString)
    for(i in 1:length(matches)){
     # i = 2
      ggg = as.numeric(unlist(matches[i])) +1
      if(levcut <= (length(ggg)-1)){
        for(j in levcut:(length(ggg)-1)){
          if(algo == "Identity"){
            ddd1 = ((Clus$seqnum[ggg[j]] - Clus$seqnum[ggg[j+1]]) / Clus$seqnum[ggg[j]]) * 100
            ddd2 = (Clus$Identity[ggg[j+1]] - Clus$Identity[ggg[j]]) * 0.5 + (Clus$Similarity[ggg[j+1]] - Clus$Similarity[ggg[j]]) * 0.5
            if(ddd1 < threshold1 || ddd2 < threshold2){
              deik = which(x$pathString == arxpath[i])
              tem = str_locate(arxpath[i],as.character(ggg[j]-1))
              fftemp = str_sub(x$pathString[deik[1]],1,tem[2])
              ff2 = str_sub(x$pathString[deik[1]],tem[1],tem[2])
              ff3 = as.numeric(ff2)
              x$clusters[deik] = ff3
              ff4 = lastlist$clep[ff3]
              ff5 = str_which(names(x),sprintf("level.%d",ff4))
              ff6 = str_which(names(x), "level.")
              if(is.na(ff5[1]) == FALSE){
                x[deik,(ff5[1]+1):max(ff6)] = NA
              }
              x$pathString[deik] = fftemp
              break()
            }
          }else if(algo == "Similarity"){
            ddd1 = ((Clus$seqnum[ggg[j]] - Clus$seqnum[ggg[j+1]]) / Clus$seqnum[ggg[j]]) * 100
            ddd2 = Clus$Similarity[ggg[j+1]] - Clus$Similarity[ggg[j]]
            #print(paste0("ddd1=",ddd1," ddd2=",ddd2))
            #print(paste0("i=",i," j=",j) )
            if(ddd1 < threshold1 || ddd2 < threshold2){
              deik = which(x$pathString == arxpath[i])
              tem = str_locate(arxpath[i],as.character(ggg[j]-1))
              fftemp = str_sub(x$pathString[deik[1]],1,tem[2])
              ff2 = str_sub(x$pathString[deik[1]],tem[1],tem[2])
              ff3 = as.numeric(ff2)
              x$clusters[deik] = ff3
              ff4 = lastlist$clep[ff3]
              ff5 = str_which(names(x),sprintf("level.%d",ff4))
              ff6 = str_which(names(x), "level.")
              if(is.na(ff5[1]) == FALSE){
                x[deik,(ff5[1]+1):max(ff6)] = NA
              }
              x$pathString[deik] = fftemp
              break()
            }
          }
        }
      }
    }
    xN <- as.Node(x)
    bo = which(lastlist$clep == lev)
    levtel = lastlist$clep
    j = 1
    jfjf = vector()
    for(i in 1:max(bo)){
      temp1 = as.numeric(names(FindNode(xN,(sprintf("%d",i)))$children))
      if(length(temp1) == 1){
        if(temp1[1] %% 2 == 1){
          ff7 = temp1[1] + 1
        }else{
          ff7 = temp1[1] - 1
        }
        ff8 = lastlist$clep[ff7]
        ff9 = str_which(names(x),sprintf("level.%d",ff8))
        deik2 = which(df[ff9] == ff7)
        x[deik2,ff9] = df[deik2,ff9]
        tem = str_locate(df$pathString[deik2[1]],as.character(ff7))
        fftemp = str_sub(df$pathString[deik2[1]],1,tem[2])
        x$pathString[deik2] = fftemp
        jfjf[j] = ff7
        j = j+1
      }
      if(is.null(FindNode(xN,(sprintf("%d",i)))$isLeaf) && is.element(i, jfjf) == FALSE){
        levtel[i] = NA
        Clus[i+1,] = NA
        ff[i+1,] = NA
      }
    }
  }
  xN <- as.Node(x)
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  plot(xN)
}

# Create custom colour scheme
cs1 = make_col_scheme(chars=c("F","W","A","I","L","V","M","C","P","G","Y","T","S","H","K","R","E","D","Q","N"),
                      cols=c("#1E90FF", "#BA55D3", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#C6E2FF", "#C6E2FF", "#FFD700", "#00EE00", "#C1FFC1", "#54FF9F", "#54FF9F", "#FF0000", "#FF0000", "#FF0000", "#FFD700", "#FFD700", "#ED9121", "#ED9121"))
#df= x
#lev = 2
#lastlist$clep = levtel
# A function to plot with logo level lev
lev = 5
seqthr = 3
idethr = 70
simthr = 
simen = FALSE
seqen = FALSE
ideen = FALSE
allen = TRUE
LogoLev <- function(lev){
  if(flagtic == TRUE) tic(sprintf("LogoLev -- Level: %d ", lev))
  if( is.element(lev,lastlist$clep) == FALSE){
    print("Den yparxei")
  }else{
    t1 = which(lastlist$clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(lastlist$clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    
    if(seqen == TRUE){
      if(length(which(Clus[t1+1,]$seqnum < seqthr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$seqnum < seqthr)]
      }
    }else if(ideen == TRUE){
      if(length(which(Clus[t1+1,]$Identity < idethr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$Identity < idethr)]
      }
    }else if (simen == TRUE){
      if(length(which(Clus[t1+1,]$Similarity < simthr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$Similarity <  simthr)]
      }
    }else if (allen == TRUE){
      if(length(which(Clus[t1+1,]$seqnum < seqthr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$seqnum < seqthr)]
      }
      if(length(which(Clus[t1+1,]$Identity < idethr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$Identity < idethr)]
      }
      if(length(which(Clus[t1+1,]$Similarity <  simthr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$Similarity <  simthr)]
      }
    }
    
    #swap(a$x,a$freq)
    #a = plyr::count(Clus[t1,]$Similarity)
    #barplot(t(as.matrix(a)), width=2, main = sprintf('Cluster.%d',cl)) 
    #legend("topright",inset=c(-0.03,0), fill=heat.colors(nrow(a)), legend=let,cex = 0.6)
    #if(length(t1) %% 3 == 0){
    #  nc = length(t1)%/%3
    #}else{
    #  nc = length(t1)%/%3 + 1
    #}
    listff = list()
    nc = 3
    for(i in 1:length(t1)){
      if(t1[i]<= max(xm2)){
        listff$temp = na.omit(df[df[names(df) == sprintf('level.%d', lastlist$clep[t1[i]])] == t1[i],]$AA.JUNCTION)
        names(listff)[length(listff)] = sprintf('Cluster.%d - num:%d - leaf', t1[i],Clus$seqnum[t1[i]+1]) 
      }else{
        x1 = as.data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", lev))] == t1[i], ]$AA.JUNCTION))
        if (nrow(x1) >0){
          names(x1)[1]= "AA.JUNCTION"
          x1 = as.character(x1$AA.JUNCTION)
          listff$temp = x1
          names(listff)[length(listff)] = sprintf('Cluster.%d - num:%d', t1[i],Clus$seqnum[t1[i]+1]) 
        } 
      }
    }
    
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    ggseqlogo(listff, ncol=nc, method = "prob",col_scheme=cs1) 
  }
}

SatLev <- function(lev){
  if(flagtic == TRUE) tic()
  if( is.element(lev,lastlist$clep) == FALSE){
    NULL
  }else{
    xN <- as.Node(df)
    t1 = which(lastlist$clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(lastlist$clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    par(mfrow=c(1,3))
    a = table(Clus[t1+1,]$seqnum)
    barplot(a, width=2, main = sprintf('Sequences for level.%d' ,lev))
    b = table(round(Clus[t1+1,]$Identity,2))
    barplot(b, width=2, main = sprintf('Identity of sequences for level.%d' ,lev))
    c = table(round(Clus[t1+1,]$Similarity,2))
    barplot(c, width=2, main = sprintf('Similarity of sequences for level.%d' ,lev))
  }
  if(flagtic == TRUE){
    en = toc(quiet = TRUE)
    cat(paste0(sprintf("LogoLev -- Level: %d ", lev),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  }
}

#df = x
# Creating a plot for cluster 5 for example
LogoCl <- function(cl){
  if(flagtic == TRUE) tic(sprintf("LogoCl -- Cluster: %d ", cl))
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  if(is.na(lastlist$clep[cl])){
    print("Den yparxei")
  }else{
    xN <- as.Node(df)
    listff = list()
    listff$temp = na.omit(df[df[names(df) == sprintf('level.%d', lastlist$clep[cl])] == cl,]$AA.JUNCTION)
    xm = as.numeric(as.data.frame(xN$leaves))
    
    if(is.element(cl,xm)){
      names(listff)= sprintf('Cluster: %s    num:%d - leaf', cl ,Clus$seqnum[cl+1])  
    }else{
      names(listff)= sprintf('Cluster: %s    num:%d', cl ,Clus$seqnum[cl+1]) 
    }
    ggseqlogo(listff, method = "prob", col_scheme=cs1)
  }
}

# A function to plot with barplot level lev
BarLev <- function(lev){
  if(flagtic == TRUE) tic(sprintf("BarLev -- Level: %d ", lev))
  if(cho == "Identity"){
    if( is.element(lev,lastlist$clep) == FALSE){
      print("Den yparxei")
    }else{
      t2 = which(lastlist$clep == lev)
      # epipleon
      xm = as.numeric(as.data.frame(xN$leaves))
      orio = min(which(lastlist$clep == lev))
      xm2 = sort(xm[xm < orio])
      t2 = sort(append(xm2,t2,after = length(xm2)))
      if(length(t2) %% 3 == 0){
        par(mfrow = c(length(t2)%/%3,3))
      }else{
        par(mfrow = c(length(t2)%/%3 + 1,3))
      }
      for(i in 1:length(t2)){
        ar = str_which(names(perlist),as.character(t2[i]))[1]
        par(xpd=TRUE)
        output <- matrix(unlist(perlist[ar]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
        if(i<= length(xm2)){
          barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d - seqnum:%d - leaf', t2[i],Clus$seqnum[t2[i]+1])) 
          
        }else{
          barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d - seqnum:%d', t2[i],Clus$seqnum[t2[i]+1])) 
        }
        legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=let,cex = 0.6)
      }
    }
  }else{
    if( is.element(lev,lastlist$clep) == FALSE){
      print("Den yparxei")
    }else{
      xN = as.Node(df)
      t2 = which(lastlist$clep == lev)
      # epipleon
      xm = as.numeric(as.data.frame(xN$leaves))
      orio = min(which(lastlist$clep == lev))
      xm2 = sort(xm[xm < orio])
      t2 = sort(append(xm2,t2,after = length(xm2)))
      if(length(t2) %% 3 == 0){
        par(mfrow = c(length(t2)%/%3,3))
      }else{
        par(mfrow = c(length(t2)%/%3 + 1,3))
      }
      for(i in 1:length(t2)){
        ar = str_which(names(persimlist),as.character(t2[i]))[1]
        par(xpd=TRUE)
        output <- matrix(unlist(persimlist[ar]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
        if(i<= length(xm2)){
          barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d - seqnum:%d - leaf', t2[i],Clus$seqnum[t2[i]+1])) 
          
        }else{
          barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d - seqnum:%d', t2[i],Clus$seqnum[t2[i]+1])) 
        }
        legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=names(sim),cex = 0.6)
      }
    }
  }
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
}
cho = "Identity"
# Create the grid (we need to miximize the pop-up window in order to have a perfect fit)
#get( getOption( "device" ) )()


# An alternative plot for cluster 3 (We must determine the cluster's permat)
BarCl <- function(cl){
  if(flagtic == TRUE) tic(sprintf("BarCl -- Cluster :%d ", cl))
  if(cho == "Identity"){
    if(is.na(lastlist$clep[cl])){
      print("Den yparxei")
    }else{
      xN <- as.Node(df)
      xm = as.numeric(as.data.frame(xN$leaves))
      par(mar=c(3,3,4,4),xpd=TRUE)
      output <- matrix(unlist(perlist[cl+1]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
      if(is.element(cl,xm)){
        barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster: %s    num:%d - leaf', cl ,Clus$seqnum[cl+1])) 
        names(listff)= sprintf('Cluster: %s    num:%d - leaf', cl ,Clus$seqnum[cl+1])  
      }else{
        barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster: %s    num:%d', cl ,Clus$seqnum[cl+1])) 
      }
      legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=let,cex = 0.6)
    }
  }else{
    if(is.na(lastlist$clep[cl])){
      print("Den yparxei")
    }else{
      xN <- as.Node(df)
      xm = as.numeric(as.data.frame(xN$leaves))
      par(mar=c(3,3,4,4),xpd=TRUE)
      output <- matrix(unlist(persimlist[cl+1]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
      if(is.element(cl,xm)){
        barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster: %s    num:%d - leaf', cl ,Clus$seqnum[cl+1])) 
        names(listff)= sprintf('Cluster: %s    num:%d - leaf', cl ,Clus$seqnum[cl+1])  
      }else{
        barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster: %s    num:%d', cl ,Clus$seqnum[cl+1])) 
      }
      legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=names(sim),cex = 0.6)
    }
  }
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
}

# Sequences and Id's for specific level
AminoLev <- function(level){
  if(flagtic == TRUE) tic(sprintf("AminoLev -- Level: %d ", level))
  if( is.element(lev,lastlist$clep) == FALSE){
    print("Den yparxei")
  }else{
    sum(Clus[Clus$level == level,]$seqnum) # akoloy8ies sto level
    x4 = data_frame("Sequence.ID" = character(0),"AA.JUNCTION" = character(0))
    gg = na.omit(unique(df[,which(names(df) == sprintf("level.%d", level))]))
    for(i in 1:length(gg)){
      clust = gg[i]
      x3 = data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$Sequence.ID), na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$AA.JUNCTION))
      names(x3) = c("Sequence.ID","AA.JUNCTION")
      se = x3$Sequence.ID
      aa = as.character(x3$AA.JUNCTION)
      aa = as.list(aa)
      names(aa) = se
      x4 = rbind(x4,x3)
    }
    se2 = x4$Sequence.ID
    aa2 = as.character(x4$AA.JUNCTION)
    aa2 = as.list(aa2)
    names(aa2) = se2
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    x4 # Print in console 
  }
}
#write.fasta(aa2, names= FALSE ,file = sprintf("Level.%d", level)) # Create a fasta file

# Sequences and Id's for specific cluster
AminoCl <- function(clust){
  #if(flagtic == TRUE) tic("AminoCl")
  #if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  if(clust == 0){
    x3 = data.frame("Sequence.ID" = df$Sequence.ID, "AA.JUNCTION" = df$AA.JUNCTION, "J.GENE.and.allele" = df$V.GENE.and.allele)
    x3
  }else{
    if(is.na(lastlist$clep[clust]) == FALSE){
     # print("Den yparxei")
    #}else{
      x3 = data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$Sequence.ID), na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$AA.JUNCTION),na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$V.GENE.and.allele))
      names(x3) = c("Sequence.ID","AA.JUNCTION","V.GENE.and.allele")
      se = x3$Sequence.ID
      aa = as.character(x3$AA.JUNCTION)
      aa = as.list(aa)
      names(aa) = se 
      x3 # Print in console 
    }
  }
}


# A function which visualize the sequences of a data frame using common letters (i.e. "A _ _ _ _ K R _ _ _ Q _ Y Y Y _ _ _ T _")
Opt <- function(df){
  if(is.list(df)){
    #print("Den yparxei")
    #if(flagtic == TRUE) tic("Opt")
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
      xar[i] = "_"
      if (f[1] == nrow(df)){
        xar[i] = names(f[1])
      }
    }
    # if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    #paste(xar,collapse = ' ')
    xar
  }
}

# Newwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
typeopt = "Identity"
OptLev <- function(lev){
  xN <- as.Node(df)
  t1 = which(lastlist$clep == lev)
  # epipleon
  xm = as.numeric(as.data.frame(xN$leaves))
  orio = min(which(lastlist$clep == lev))
  xm2 = sort(xm[xm < orio])
  t1 = sort(append(xm2,t1,after = length(xm2)))
  if(typeopt == "Identity"){
    mm = as.data.frame(matrix(0,nrow=length(t1), ncol=str_length(udata$AA.JUNCTION[1])) )
    for (i in 1:length(t1)) {
      mm[i,] = Opt(AminoCl(t1[i]))
      if(i<= length(xm2)){
        row.names(mm)[i] = sprintf('cluster.%d - seqnum:%d - leaf',t1[i],Clus$seqnum[t1[i]+1])
      }else{
        row.names(mm)[i] = sprintf('cluster.%d - seqnum:%d ',t1[i],Clus$seqnum[t1[i]+1])
      }
    }
  }else if (typeopt == "IdentityFreq"){
    mm = as.data.frame(matrix(0,nrow=2*length(t1), ncol=str_length(udata$AA.JUNCTION[1])) )
    temp = 0
    for (i in 1:length(t1)) {
      temp = temp + 1
      mm[temp,] = OptNew(AminoCl(t1[i]))[1,]
      if(i<= length(xm2)){
        row.names(mm)[temp] = sprintf('cluster.%d - seqnum:%d - leaf',t1[i],Clus$seqnum[t1[i]+1])
      }else{
        row.names(mm)[temp] = sprintf('cluster.%d - seqnum:%d ',t1[i],Clus$seqnum[t1[i]+1])
      }
      temp = temp +1
      mm[temp,] = OptNew(AminoCl(t1[i]))[2,]
      row.names(mm)[temp] = sprintf('cluster.%d - percentage',t1[i])
    }
  }else if (typeopt == "Similarity"){
    mm = as.data.frame(matrix(0,nrow=length(t1), ncol=str_length(udata$AA.JUNCTION[1])) )
    for (i in 1:length(t1)) {
      mm[i,] = Opt2(AminoCl(t1[i]),altsim)
      if(i<= length(xm2)){
        row.names(mm)[i] = sprintf('cluster.%d - seqnum:%d - leaf',t1[i],Clus$seqnum[t1[i]+1])
      }else{
        row.names(mm)[i] = sprintf('cluster.%d - seqnum:%d ',t1[i],Clus$seqnum[t1[i]+1])
      }
    }
  }else if(typeopt == "SimilarityFreq"){
    mm = as.data.frame(matrix(0,nrow=2*length(t1), ncol=str_length(udata$AA.JUNCTION[1])) )
    temp = 0
    for (i in 1:length(t1)) {
      temp = temp + 1
      mm[temp,] = Opt2New(AminoCl(t1[1]))[1,]
      if(i<= length(xm2)){
        row.names(mm)[temp] = sprintf('cluster.%d - seqnum:%d - leaf',t1[i],Clus$seqnum[t1[i]+1])
      }else{
        row.names(mm)[temp] = sprintf('cluster.%d - seqnum:%d ',t1[i],Clus$seqnum[t1[i]+1])
      }
      temp = temp +1
      mm[temp,] = Opt2New(AminoCl(t1[1]))[2,]
      row.names(mm)[temp] = sprintf('cluster.%d - percentage',t1[i])
    }
  }
  #fgh
  mm
}





#Opt(AminoCl(1))
#fgh = sprintf('cluster.%d:    %s',t1[1],(Opt(AminoCl(t1[1]))))
#for(i in 2:length(t1)){
#  #fgh = Opt(AminoCl(t1[i]))
#  fgh = paste(fgh,sprintf('cluster.%d:    %s',t1[i],(Opt(AminoCl(t1[i])))),sep = "\n")
#}
#cncn = cat(fgh,sep = "\n")
poso = TRUE
df = AminoCl(17)
OptNew <- function(df){
  if(is.list(df)){
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    xarpos <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
      if (length(which(f == max(f))) > 1){
        tempopt = "-"
        tempopt = paste(tempopt,names(which(f == max(f)))[1],sep = "")
        #tempopt =  names(which(f == max(f)))[1]
        for(j in 2:(length(which(f == max(f))))){
          tempopt = paste(tempopt, names(which(f == max(f)))[j],sep = "|")
        }
        tempopt = paste(tempopt,"-",sep = "")
        xar[i] = tempopt
        xarpos[i] = round(((max(f) / nrow(df)) * 100), digits = 2)
      }else{
        xar[i] = names(which(f == max(f)))
        xarpos[i] = round(((max(f) / nrow(df)) * 100), digits = 2)
      }
    }
    #max1 = max(str_length(xar))
    #max2 = max(str_length(xarpos))
    #max3 = max(max1,max2)
    #if(poso == TRUE){
    #  paste(str_pad(xarpos,max3),collapse = ' ')
    #}else{
    #  paste(str_pad(xar,max3),collapse = ' ')
    #}
    rbind(xar,xarpos)
  }
  
}
#OptNew(AminoCl(1))
mx = as.data.frame(rbind(xar,xarpos))
xmmm = t(rbind(xar,xarpos))

# A function which visualize the sequences of a data frame using similarity groups (i.e. "A _ _ _ _ K Am _ Ba _ Ac _ Y Y Y _ _ _ T _")
Opt2 <- function(df,alt){
  if(is.list(df)){
    #if(flagtic == TRUE) tic("Opt2")
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
          if(alt == TRUE){
            xar[i] = names(altsim[d])
          }else{
            xar[i] = names(sim[d])
          }
        } 
      }
      
    }
    #if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    paste(xar,collapse = ' ') 
  }
}
#Opt2(AminoCl(5),TRUE)

Opt2New <- function(df){
  if(is.list(df)){
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    xarpos <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
      for(j in 1:length(f)){
        if( length(which(str_detect(names(f),names(sim[str_which(sim,names(f)[j])])))) == 0 ){
          names(f)[j] = names(sim[str_which(sim,names(f)[j])])
        }else if(which(str_detect(names(f),names(sim[str_which(sim,names(f)[j])]))) != j){
          f[which(str_detect(names(f),names(sim[str_which(sim,names(f)[j])])))] = f[which(str_detect(names(f),names(sim[str_which(sim,names(f)[j])])))] + f[j]
          names(f)[j] = "dip"
        }
      }
      if(length(which(names(f) == "dip")) != 0 ){
        f = f[-which(names(f) == "dip")] 
      }
      if (length(which(f == max(f))) > 1){
        tempopt = "-"
        tempopt = paste(tempopt,names(which(f == max(f)))[1],sep = "")
        for(j in 2:(length(which(f == max(f))))){
          tempopt = paste(tempopt, names(which(f == max(f)))[j],sep = "|")
        }
        tempopt = paste(tempopt,"-",sep = "")
        xar[i] = tempopt
        xarpos[i] = round(((max(f) / nrow(df)) * 100), digits = 2)
      }else{
        xar[i] = names(which(f == max(f)))
        xarpos[i] = round(((max(f) / nrow(df)) * 100), digits = 2)
      }
    }
    #paste(xar,collapse = ' ')
    #if(poso == TRUE){
    #  paste(xarpos,collapse = ' ')
    #}else{
    #  paste(xar,collapse = ' ')
    #}
    rbind(xar,xarpos)
  }
}
#Opt2New(AminoCl(5))

Id <- function(ff){
  if(flagtic == TRUE) tic(sprintf("Id -- Type: $s ", cho))
  if(cho == "Identity"){
    par(xpd=TRUE)
    pp <- data.frame(x = ff$level,y = ff$sumper)
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    ggplot(na.omit(pp),aes(x = x,y = y)) + stat_sum()
  }else{
    par(xpd=TRUE)
    pp <- data.frame(x = ff$level,y = ff$sumper2)
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    ggplot(na.omit(pp),aes(x = x,y = y)) + stat_sum()
  }
}


cc <- function(df){
  if(flagtic == TRUE) tic("Collapsible Tree")
  df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
  ii = Clus$seqnum
  ii = as.character(ii)
  
  trid1 = df[str_which(names(df), "level.")] # Data frame with identities percentage
  kk1 = trid1
  trsim = df[str_which(names(df), "level.")] # Data frame with similarities percentage
  iii = df[str_which(names(df), "level.")]
  for(i in 1:length(trid1)){
    tempg = trid1[,i]+1
    trid1[i] = round(Clus$Identity[tempg],digits = 2)
    temph = trsim[,i]+1
    trsim[i] = round(Clus$Similarity[temph], digits = 2)
    tempi = iii[,i] + 1
    iii[i] = ii[tempi]
  }
  
  for(i in 2:length(str_which(names(df), "level."))){
    trid1[,i] = paste(kk1[,i],iii[,i],trid1[,i], trsim[,i],sep = " ")
    for(j in 1:nrow(df)){
      if(str_detect(trid1[j,i],"NA") == TRUE){
        trid1[j,i] = NA
      }
    }
  }
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  collapsibleTree(
    trid1,
    hierarchy = colnames(df[str_which(names(df), "level.")]),
    fill = c("jj", na.omit(ii)),
    width = 1820,
    height = 775,
    collapsed = FALSE
  )
}

idenlev <- function(lev){
  if(flagtic == TRUE) tic(sprintf("idenlev -- Level: %d ", lev))
  if(cho == "Identity"){
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    t1 = which(lastlist$clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(lastlist$clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    sum(na.omit(Clus[t1+1,]$Identity)) /  length(na.omit(Clus[t1+1,]$Identity))
  }else{
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    t1 = which(lastlist$clep == lev)
    # epipleon
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(lastlist$clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    sum(na.omit(Clus[t1+1,]$Similarity)) /  length(na.omit(Clus[t1+1,]$Similarity))
  }
}

idencl <- function(cl){
  if(flagtic == TRUE) tic(sprintf("idencl -- Cluster :%d ", cl))
  if(is.element(cl,Clus$ClusterId)){
    if(cho == "Identity"){
      if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
      na.omit(Clus[Clus$ClusterId == cl,])$Identity
    }else{
      if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
      na.omit(Clus[Clus$ClusterId == cl,])$Similarity 
    }
  }else{
    "Den yparxei"
  } 
}

EmPin <- function(lastlist,Clus,dfsd,flagtic){
  if(flagtic == TRUE) tic("EmPin")
  
  for (i in 0:(max(na.omit(lastlist$clep)))) {
    if(i == 0){
      t1 = Clus$Identity[i+1]
      t2 = sd(Clus$Identity[i+1])
      t3 = Clus$Similarity[i+1]
      t4 = sd(Clus$Similarity)[i+1]
    }else{
      t = which(lastlist$clep == i)
      # epipleon
      xm = as.numeric(as.data.frame(xN$leaves))
      orio = min(which(lastlist$clep == i))
      xm2 = sort(xm[xm < orio])
      t = sort(append(xm2,t,after = length(xm2)))
      t1 = sum(na.omit(Clus[t+1,]$Identity)) /  length(na.omit(Clus[t+1,]$Identity))
      t2 = sd(na.omit(Clus[t+1,]$Identity))
      t3 = sum(na.omit(Clus[t+1,]$Similarity)) /  length(na.omit(Clus[t+1,]$Similarity))
      t4 = sd(na.omit(Clus[t+1,]$Similarity))
    }
    rows = sprintf("level.%d", i)
    dfsd[i+1,] = c(rows,t1,t2,t3,t4)
  }
  colnames(dfsd)[1]="Level"
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  dfsd
}

lastlist = lastlist2
df = lastlist$udata
lev = 5
thrt = "Distance"
thr = 0
netyp = "Whole"
enthr = TRUE
threshold1 = 7
threshold2 = 10
net_sil = FALSE
Clus = Clus2
Netw <- function(lev,thr,thrt,netyp,df,lastlist,Clus,threshold1,threshold2,enthr){
  if(flagtic == TRUE) tic(sprintf("Network -- level: %d, Network type: %s, Threshold type: %s, Threshold value: %d, Using silhouette: %d ", lev,netyp,thrt,thr,net_sil))
  df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
  df_args <- c(df[str_which(names(df), "level.")], sep="/")
  if(lev == max(na.omit(lastlist$clep))){
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
  
  if(enthr == TRUE){
    matches <- regmatches(unique(x$pathString), gregexpr("[[:digit:]]+", unique(x$pathString)))
    arxpath = unique(x$pathString)
    for(i in 1:length(matches)){
      ggg = as.numeric(unlist(matches[i])) +1
      if(levcut <= (length(ggg)-1)){
        for(j in levcut:(length(ggg)-1)){
          if(algo == "Identity"){
            ddd1 = ((Clus$seqnum[ggg[j]] - Clus$seqnum[ggg[j+1]]) / Clus$seqnum[ggg[j]]) * 100
            ddd2 = (Clus$Identity[ggg[j+1]] - Clus$Identity[ggg[j]]) * 0.5 + (Clus$Similarity[ggg[j+1]] - Clus$Similarity[ggg[j]]) * 0.5
            if(ddd1 < threshold1 || ddd2 < threshold2){
              deik = which(x$pathString == arxpath[i])
              tem = str_locate(arxpath[i],as.character(ggg[j]-1))
              fftemp = str_sub(x$pathString[deik[1]],1,tem[2])
              ff2 = str_sub(x$pathString[deik[1]],tem[1],tem[2])
              ff3 = as.numeric(ff2)
              x$clusters[deik] = ff3
              ff4 = lastlist$clep[ff3]
              ff5 = str_which(names(x),sprintf("level.%d",ff4))
              ff6 = str_which(names(x), "level.")
              if(is.na(ff5[1]) == FALSE){
                x[deik,(ff5[1]+1):max(ff6)] = NA
              }
              x$pathString[deik] = fftemp
              break()
            }
          }else if(algo == "Similarity"){
            ddd1 = ((Clus$seqnum[ggg[j]] - Clus$seqnum[ggg[j+1]]) / Clus$seqnum[ggg[j]]) * 100
            ddd2 = Clus$Similarity[ggg[j+1]] - Clus$Similarity[ggg[j]]
            print(paste0("ddd1=",ddd1," ddd2=",ddd2))
            print(paste0("i=",i," j=",j) )
            if(ddd1 < threshold1 || ddd2 < threshold2){
              deik = which(x$pathString == arxpath[i])
              tem = str_locate(arxpath[i],as.character(ggg[j]-1))
              fftemp = str_sub(x$pathString[deik[1]],1,tem[2])
              ff2 = str_sub(x$pathString[deik[1]],tem[1],tem[2])
              ff3 = as.numeric(ff2)
              x$clusters[deik] = ff3
              ff4 = lastlist$clep[ff3]
              ff5 = str_which(names(x),sprintf("level.%d",ff4))
              ff6 = str_which(names(x), "level.")
              if(is.na(ff5[1]) == FALSE){
                x[deik,(ff5[1]+1):max(ff6)] = NA
              }
              x$pathString[deik] = fftemp
              break()
            }
          }
        }
      }
    }
    xN <- as.Node(x)
    bo = which(lastlist$clep == lev)
    levtel = lastlist$clep
    j = 1
    jfjf = vector()
    for(i in 1:max(bo)){
      temp1 = as.numeric(names(FindNode(xN,(sprintf("%d",i)))$children))
      if(length(temp1) == 1){
        if(temp1[1] %% 2 == 1){
          ff7 = temp1[1] + 1
        }else{
          ff7 = temp1[1] - 1
        }
        ff8 = lastlist$clep[ff7]
        ff9 = str_which(names(x),sprintf("level.%d",ff8))
        deik2 = which(df[ff9] == ff7)
        x[deik2,ff9] = df[deik2,ff9]
        tem = str_locate(df$pathString[deik2[1]],as.character(ff7))
        fftemp = str_sub(df$pathString[deik2[1]],1,tem[2])
        x$pathString[deik2] = fftemp
        jfjf[j] = ff7
        j = j+1
      }
      if(is.null(FindNode(xN,(sprintf("%d",i)))$isLeaf) && is.element(i, jfjf) == FALSE){
        levtel[i] = NA
        Clus[i+1,] = NA
        ff[i+1,] = NA
      }
    }
    xN <- as.Node(x)
    
    #df = x
    #lastlist$udata = df
    #lastlist$clep = levtel
    
    sss <- regmatches(unique(x$pathString), gregexpr("[[:digit:]]+", unique(x$pathString)))
    arxpath = unique(x$pathString)
    for(i in 1:length(sss)){
      sss1 = as.numeric(unlist(sss[i]))
      if (length(sss1) > 3){
        deik = which(x$pathString == arxpath[i])
        tem = str_locate(arxpath[i],as.character(rev(sss1)[3]))
        fftemp = str_sub(arxpath[i],tem[1],str_length(arxpath[i]))
        x$pathString[deik] = fftemp 
      }
    }
  }
  
  
  #xN <- as.Node(x)
  
  bo = which(lastlist$clep == lev)
  ffg = as.data.frame(matrix(0,nrow = max(bo)+1,ncol = (max(bo)+1))) # the distance between 2 Nodes
  ffg2 = as.data.frame(matrix(0,nrow = max(bo)+1,ncol = (max(bo)+1)))
  ffg2sim = as.data.frame(matrix(0,nrow = max(bo)+1,ncol = (max(bo)+1)))
  ffg3 = as.data.frame(matrix(NA,nrow = max(bo)+1,ncol = (max(bo)+1))) # the final distance combining ffg and ffg2
  
  
  if(enthr == TRUE){
    bfbf = unique(as.numeric(unlist(sss)))
    
    #tic()
    s = sapply(1:length(bfbf), function(i){
      hhh = str_which(sss,sprintf("%d", bfbf[i]))
      tempN1 = unlist(sss[hhh[1]])
      tempN1 = tempN1[1:str_which(tempN1,sprintf("%d", bfbf[i]))]
      temp2N1 = paste(Opt(AminoCl(bfbf[i])),collapse = ' ')
      temp3N1 = paste(Opt2(AminoCl(bfbf[i]),TRUE),collapse = ' ')
      sapply(i:length(bfbf),function(j) {
        hhh2 = str_which(sss,sprintf("%d", bfbf[j]))
        tempN2 = unlist(sss[hhh2[1]])
        tempN2 = tempN2[1:str_which(tempN2,sprintf("%d", bfbf[j]))]
        tem = which(tempN1 == tempN2)
        d1 = length(tempN1) - tem[length(tem)]
        d2 = length(tempN2) - tem[length(tem)]
        temp2N2 = paste(Opt(AminoCl( bfbf[j])),collapse = ' ')
        temp2N2 = str_replace_all(temp2N2,"_"," ")
        temp3N2 = paste(Opt2(AminoCl( bfbf[j]),TRUE),collapse = ' ')
        temp3N2 = str_replace_all(temp3N2,"_"," ")
        set(ffg,  bfbf[i]+1,  bfbf[j]+1, d1 + d2) 
        set(ffg2, (bfbf[i]+1),  bfbf[j]+1, stringdist(temp2N1,temp2N2,method = "lv" )) 
        set(ffg2sim, (bfbf[i]+1),  bfbf[j]+1, stringdist(temp3N1,temp3N2,method = "lv" )) 
      })
    })
    #toc()
  }else{
    #tic()
    s = sapply(1:(max(bo)+1), function(i){
      tempN1 = FindNode(xN,(sprintf("%d", i-1)))
      tempN1 = tempN1$path
      temp2N1 = paste(Opt(AminoCl(i-1)),collapse = ' ')
      temp3N1 =  paste(Opt2(AminoCl(i-1),TRUE),collapse = ' ')
      sapply(i:(max(bo)+1),function(j) {
        tempN2 = FindNode(xN,(sprintf("%d", j-1)))
        tempN2 = tempN2$path
        tem = which(tempN1 == tempN2)
        d1 = length(tempN1) - tem[length(tem)]
        d2 = length(tempN2) - tem[length(tem)]
        temp2N2 = paste(Opt(AminoCl(j-1)),collapse = ' ')
        temp2N2 = str_replace_all(temp2N2,"_"," ")
        temp3N2 = paste(Opt2(AminoCl(j-1),TRUE),collapse = ' ')
        temp3N2 = str_replace_all(temp3N2,"_"," ")
        set(ffg, i, j, as.integer(d1 + d2)) 
        set(ffg2, i, j, stringdist(temp2N1,temp2N2,method = "lv" )) 
        set(ffg2sim, i, j, stringdist(temp3N1,temp3N2,method = "lv" )) 
      })
    })
    #toc()
  }
  
  
  
  #i = 1
  # anw trigwnikos
  
  
  #find max
  ma = max(ffg[is.na(ffg) == FALSE])
  ffg[ffg == 0] = ma
  ffg2[ffg2 == 0] = str_length(lastlist$udata$AA.JUNCTION[1])
  ffg2sim[ffg2sim == 0] = str_length(lastlist$udata$AA.JUNCTION[1])
  ffg[lower.tri(ffg)] = NA
  ffg2[lower.tri(ffg2)] = NA
  ffg2sim[lower.tri(ffg2sim)] = NA
  #ffg[ffg == 0] = NA
  #ffg2[ffg2 == 0] = NA
  #ffg2sim[ffg2sim == 0] = NA
  #ffg[which(ffg == 0)] = ma #gia na mhn fainontai velakia pisw
  ffg2 = str_length(lastlist$udata$AA.JUNCTION[1]) - ffg2
  ffg2sim = str_length(lastlist$udata$AA.JUNCTION[1]) - ffg2sim
  ffg3 = (2/str_length(lastlist$udata$AA.JUNCTION[1])) * ffg2 + (2/str_length(lastlist$udata$AA.JUNCTION[1])) * ffg2sim + (2 / ma) * (ma - ffg)
  
  diag(ffg) = 0
  diag(ffg2) = 0
  diag(ffg2sim) = 0
  diag(ffg3) = 0
  
  normalize <- function(x) {
    return ((x - min(x[is.na(x) == FALSE])) / (max(x[is.na(x) == FALSE]) - min(x[is.na(x) == FALSE])))
  }
  
  if(thrt == "Distance"){
    thrtyp = ffg
  }else if (thrt == "StrIdentSimilarity"){
    thrtyp = ffg2
  }else if (thrt == "StrGroupSimilarity"){
    thrtyp = ffg2sim
  }else{
    thrtyp = ffg3
  }
  
  tempor = ffg
  tempor2 = ffg2
  tempor2sim = ffg2sim
  tempor3 = ffg3
  tempor[thrtyp < thr] = 0
  tempor2[thrtyp > thr] = 0
  tempor2sim[thrtyp > thr] = 0
  tempor3[thrtyp> thr] = 0
  
  if(thrt == "Distance"){
    thrtyp2 = tempor
  }else if (thrt == "StrIdentSimilarity"){
    thrtyp2 = tempor2
  }else if (thrt == "StrGroupSimilarity"){
    thrtyp2 = tempor2sim
  }else{
    thrtyp2 = tempor3
  }
  
  if(net_sil == TRUE){
    newffg = normalize(ffg)
    newffg2 = normalize(ffg2)
    newffg2sim = normalize(ffg2sim)
    newffg3 = newffg * 0.5 + newffg2 * 0.25 + newffg2sim * 0.25  
    nnnn = normalize(newffg3)
    
    jhj = silhouette(Clus$level[1:(max(bo)+1)],t(nnnn))
    matches <- regmatches(unique(x$pathString), gregexpr("[[:digit:]]+", unique(x$pathString)))
    tttsyn = lapply(1:length(matches),function(i){
      ttt = jhj[as.numeric(unlist(matches[i])),3]
      sapply(length(ttt):1,function(j){
        if(j > 1){
          if(ttt[j] >= ttt[j-1]){
            set(tempor3, j-1, j, 0)
            set(ffg3, j-1, j, 0)
            set(thrtyp, j-1, j, 0)
            set(thrtyp2, j-1, j, 0)
          }
        }
      })
    })
  }
  
  jjj = matrix(1,nrow = max(bo)+1,ncol = (max(bo)+1))
  #jjj = matrix(1,nrow = length(bb),ncol = length(bb))
  net0 <- graph_from_adjacency_matrix(jjj,mode = "upper")
  net1 <- graph_from_adjacency_matrix(jjj,mode = "upper")
  
  bb = vector(length = (max(bo)+1))
  bb[1]=0
  bb[2:length(bb)] = lastlist$clep[1:max(bo)]
  #bb = na.omit(bb) 
  #XRWMA
  # Generate colors based on media type:
  colrs <- c("#1E90FF", "#BA55D3", "#0000FF", "#557fd2", "#54d17e", "#8aad62", "#C6E2FF", "#e5e234", "#FFD700", "#00EE00", "#C1FFC1", "#ea8509", "#54FF9F", "#FF0000", "#ed3b1c", "#ed1c7a", "#0c0c0c", "#b8d8af", "#ED9121","#45f713")
  colrs = colrs[1:(max(na.omit(bb) )+1)]
  V(net0)$color <- colrs[bb+1]
  V(net1)$color <- colrs[bb+1]
  
  # Compute node degrees (#links) and use that to set node size:
  deg <- log(Clus$seqnum[1:length(bb)] / Clus$seqnum[1] * 100) * 10
  V(net0)$size <- deg
  V(net1)$size <- deg
  
  # The labels are currently node IDs.
  # Setting them to NA will render no labels:
  V(net0)$label <- Clus$ClusterId[1:length(bb)]
  V(net1)$label <- Clus$ClusterId[1:length(bb)]
  
  # Set edge width based on weight:
  hhh = na.omit(as.vector(t(ffg3)))
  hhh1 = na.omit(as.vector(t(tempor3)))
  E(net0)$width <- hhh 
  E(net1)$width <- hhh1 
  
  #change arrow size and edge color:
  E(net0)$arrow.size <- .2
  E(net0)$edge.color <- "gray80"
  E(net1)$arrow.size <- .2
  E(net1)$edge.color <- "gray80"
  pal1 <- rainbow(6, alpha=1) 
  
  net0.copy <- igraph::delete.edges(net0, which(E(net0)$width == 0))
  #net0.copy <- igraph::delete.edges(net0, is.na(V(net0)$size))
  net0.copy <- igraph::delete.vertices(net0.copy, is.na(V(net0)$size))
  net1.copy <- igraph::delete.edges(net1, which(E(net1)$width == 0))
  
  plot(net0.copy, edge.color=na.omit(pal1[( E(net0.copy)$width %/% 1) +1]), edge.curved=.1, vertex.label.color = "black") #plot the network graph
  legend("topleft", inset=c(0.1,0.2), colnames(df[str_which(names(df), "level.")])[1:(max(na.omit(bb))+1)], pch=21,
         col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
  legend("topright", inset=c(0.1,0.2), c("0-1","1-2","2-3","3-4","4-5","5"), pch=21,
         col="#777777", pt.bg=pal1, pt.cex=2, cex=.8, bty="n", ncol=1,title = "Relationship strength")
  matrix.heatmap(thrtyp)
  
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  listnet = list("net0.copy" = net0.copy, "net1.copy" = net1.copy, "bb" = bb, "colrs" = colrs, "ma" = ma, "ffg" = ffg, "ffg2" = ffg2, "ffg2sim" = ffg2sim, "ffg3" = ffg3 ,"tempor" = tempor, "tempor2" = tempor2, "tempor2sim" = tempor2sim, "tempor3" = tempor3,"thrtyp" = thrtyp, "thrtyp2" = thrtyp2)
  return(listnet)
}

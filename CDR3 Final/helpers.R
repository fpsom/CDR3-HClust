######################################## Algorithm Implementation ########################################
# A Multi-metric Algorithm for Hierarchical Clustering of Same-length Protein Sequences - Shiny helper.R
# Tsarouchis Sotrios - Filippos
# email: sotirisftsar@gmail.com
# AEM: 7999

######################################## Function Matrices ########################################
# Function Matrices compute the apsolute matrix, the matrix with percentages, the absolute similarity matrix and the percentage similarity matrix
# Function Matrices contains the old function Finish. Finish compute cluster's identity and similarity and stop cluster's division if identity or similarity is greater than endper
Matrices <- function(list1,leaf,let,sim,d,algo,algocol,backcol,backcolj6,logFile){
  # Count the execution time of Matrices
  tic()
  # Keep data up to date
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
  
  # Find the sequences with gene J6
  ind1 = str_which(udata[udata$clusters == br,]$J.GENE.and.allele,"J6")
  qw = 1:length(udata[udata$clusters == br,]$AA.JUNCTION)
  if(length(ind1) == 0){
    ind2 = qw
  }else{
    ind2 = qw[-ind1]
  }
  
  # Initialize matrices for J6 sequences
  mymat1 = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
  simmat1 = matrix(0,nrow = length(sim),ncol =str_length(udata$AA.JUNCTION[1]))
  rownames(mymat1) = let
  rownames(simmat1) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
  
  # Initialize matrices for NO J6 sequences
  mymat2 = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
  simmat2 = matrix(0,nrow = length(sim),ncol =str_length(udata$AA.JUNCTION[1]))
  rownames(mymat2) = let
  rownames(simmat2) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
  
  # Initialize matrices for all sequences
  permat = matrix(0,nrow=length(let) + 1, ncol=str_length(udata$AA.JUNCTION[1]))
  rownames(permat) = c(let,"Entropy")
  persim = matrix(0,nrow = length(sim)+1,ncol =str_length(udata$AA.JUNCTION[1]))
  rownames(persim) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am","Entropy")

  # Compute matrices for J6 sequences
  if(length(ind1) != 0){
    trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION[ind1],"")
    align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
    for(i in (algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)){
      temptab = plyr::count(align[i],vars = colnames(align)[i])
      names(temptab)[1] = "X1"
      mymat1[which(is.na(match(names(mymat1[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
      temptab1 = temptab
      gg = as.vector(unlist(temptab1[1]))
      temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),gg[x]))])))
      temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
      simmat1[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim) == x)),i] = as.vector(unlist(temptab1[2]))
    }
  }
  
  # Compute matrices for NO J6 sequences
  if(length(ind2) != 0 ){
    trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION[ind2],"")
    align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
    for(i in (algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)){
      temptab = plyr::count(align[i],vars = colnames(align)[i])
      names(temptab)[1] = "X1"
      mymat2[which(is.na(match(names(mymat2[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
      permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
      temptab1 = temptab
      gg = as.vector(unlist(temptab1[1]))
      temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),gg[x]))])))
      temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
      simmat2[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim) == x)),i] = as.vector(unlist(temptab1[2]))
      persim[length(sim)+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
    }
  }
  
  # Combine matrices for No J6 and J6 sequences
  mymat = mymat1 + mymat2
  simmat = simmat1 + simmat2
  meg = length(ind1) + length(ind2)
  if(length(ind2) == 0 || backcol == backcolj6){
    permat[1:length(let),(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (mymat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
    persim[1:length(sim),(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (simmat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
  }else{
    permat[1:length(let),(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (mymat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
    persim[1:length(sim),(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] = (simmat[,(algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
    if(backcolj6 != 0 ){
      permat[1:length(let),(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] = (mymat[,(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] / length(ind2)) * 100
      persim[1:length(sim),(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] = (simmat[,(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcol)] / length(ind2)) * 100
    }
  }
  
  # Find the Entropy for all the sequences
  trimudata = strsplit(udata[udata$clusters == br,]$AA.JUNCTION,"")
  align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
  for(i in (algocol+1):(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]) - backcolj6)){
    temptab = plyr::count(align[i],vars = colnames(align)[i])
    names(temptab)[1] = "X1"
    permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
    temptab1 = temptab
    gg = as.vector(unlist(temptab1[1]))
    temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),gg[x]))])))
    temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
    persim[length(sim)+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
  }
  
  if(algocol != 0 && backcol!=0){
    # Initialize matrices
    mymat3 = matrix(0,nrow=length(let), ncol=str_length(udata$AA.JUNCTION[1]))
    simmat3 = matrix(0,nrow = length(sim),ncol =str_length(udata$AA.JUNCTION[1]))
    rownames(mymat3) = let
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
      gg = as.vector(unlist(temptab1[1]))
      temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),gg[x]))])))
      temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
      simmat3[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim) == x)),exc[i]] = as.vector(unlist(temptab1[2]))
      persim3[length(sim)+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
    }
    permat3[1:length(let),] = (mymat3 / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100
    persim3[1:length(sim),] = (simmat3 / length(udata[udata$clusters == br,]$AA.JUNCTION)) * 100
    permat[1:length(let),exc] = permat3[1:length(let),exc] 
    persim[1:length(sim),exc] = persim3[1:length(sim),exc] 
  }
  
  # Keep a list with total permat matrix
  listxx$temp = permat
  names(listxx)[length(listxx)] = sprintf('permat_br.%d', br) # Save the permat with this format
  # Keep a list with total persim matrix
  listyy$temp = persim
  names(listyy)[length(listyy)] = sprintf('persim_br.%d', br) # Save the persim with this format
  
  ################################ Old Finish ###########################
  t1 = which(permat[-nrow(permat),] == 100,arr.ind = TRUE)
  sumper = (length(as.numeric(t1[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
  t2 = which(persim[-nrow(persim),] == 100,arr.ind = TRUE)
  sumper2 = (length(as.numeric(t2[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
  vv = length(udata[udata$clusters == br,]$AA.JUNCTION)
  dfsum[nrow(dfsum) + 1,] = c(sumper,sumper2,br,vv)
  
  # Check if the clusters division needs to stop
  if(algo == "Identity"){
    if (sumper > endper){
      nn = TRUE  # When nn = TRUE the percentage of sumper < endper%
      if (sumper2 > endper){
        leaf2 = TRUE
      }
    }
  }else{
    if (sumper2 > endper){
      leaf2 = TRUE
    }
  }
  
  # The return list
  result1 = list("listax" = list1$listax,"ggdf" = list1$ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep"= list1$ep, "clep" = list1$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "ela" = list1$ela, "cel" = list1$cel,"endper" = list1$endper, "last" = last, "progend" = progend, "leaf" = leaf, "leaf2" = leaf2)
  en = toc(quiet = TRUE)
  cat(paste0("Matrices","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  if ( leaf == TRUE){
    nn = FALSE
    result1 = list("listax" = list1$listax,"ggdf" = list1$ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep"= list1$ep, "clep" = list1$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "ela" = list1$ela, "cel" = list1$cel,"endper" = list1$endper, "last" = last, "progend" = progend, "leaf" = leaf, "leaf2" = leaf2)
    return(result1)
  }else{
    return(result1)
  }
}


######################################## Function Choice ########################################
# Function Choice choose which matrix cell will be used for the division of the data 
# Function Choice contains the old Divide. Divide is responsible for dividing the sequences based on the cell given.
# Function Choice also contains the old Control. Control is responsible for updating level, cluster and other counters as well as for stopping the algorithm 
Choice <- function(list2,leaf,let,sim,d,algo,algocol,backcol,backcolj6,logFile){
  # Count the execution time of Choices
  tic()
  # Keep data up to date
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
  
  #Find desired cell for division
  if(nn == TRUE || (algo == "Similarity") ){
    permat = list2$persim # If sumper < endper% we want to check only the persim matrix
  }else{
    permat =list2$permat  # Else permat and if it is necessary the persim matrix
    persim = list2$persim
  }
  cel = which(permat[,(algocol + 1):(ncol(permat)-backcol)] == max(permat[,(algocol + 1):(ncol(permat)-backcol)]), arr.ind = TRUE)
  ela = 1
  poss = max(permat[,(algocol + 1):(ncol(permat)-backcol)])
  # We exclude the 100 % from the max values
  if (max(permat[,(algocol + 1):(ncol(permat)-backcol)]) == 100){ 
    cel = which(permat[,(algocol + 1):(ncol(permat)-backcol)] == max(permat[,(algocol + 1):(ncol(permat)-backcol)][permat[,(algocol + 1):(ncol(permat)-backcol)]!=max(permat[,(algocol + 1):(ncol(permat)-backcol)])]), arr.ind = TRUE) # The desired cell
    poss = max(permat[,(algocol + 1):(ncol(permat)-backcol)][permat[,(algocol + 1):(ncol(permat)-backcol)]!=max(permat[,(algocol + 1):(ncol(permat)-backcol)])])
  }
  
  # Check if clusters division is stopped or not
  if(leaf2 == FALSE && poss != 0){
    # if cel contain more than one cells, find the best cell matching some criteria
    ela = 1
    if(algo == "Identity"){
      if ((length(cel)/2) > 1){
        dddff = min(permat[,(algocol + 1):(ncol(permat)-backcol)][nrow(permat[,(algocol + 1):(ncol(permat)-backcol)]),cel[,2]])
        dddff2 = which(permat[,(algocol + 1):(ncol(permat)-backcol)][nrow(permat[,(algocol + 1):(ncol(permat)-backcol)]),cel[,2]] == dddff)
        ela = dddff2[1]
        if (length(dddff2)>1 && nn == FALSE){ # If the vector has 2 or more numbers means that we have columns with the same entropy and nn = FALSE in order not to double check the persim
          dddff3 = persim[,(algocol + 1):(ncol(permat)-backcol)][sapply(1:length(dddff2), function (x){ str_which(names(sim),gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),let[cel[x,1]]))])))}),cel[,2]]
          dddff3 = max(diag(dddff3))
          dddff4 = which(diag( persim[,(algocol + 1):(ncol(permat)-backcol)][sapply(1:length(dddff2), function (x){str_which(names(sim),gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),let[cel[x,1]]))])))}),cel[,2]]) == dddff3)
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
    
    ##################################################### Old Divide ############################
    # If we need a new level, then we create a new column in udata dataframe with its name (level.ep)
    if (met == 0){
      udata$temp = NA
      names(udata)[length(udata)] = sprintf('level.%d', ep)
    }
    
    if (algo == "Identity"){
      # Find sequences contains the cell's letter in cell's position 
      x1 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), let[cel[ela,1]]), "TRUE")
      mk1 = length(x1)
      cltp1 = cl + 1
      y1 = udata[udata$clusters == br,]$AA.JUNCTION
      z1 = y1[x1]
      
      # The other sequences of the cluster
      x2 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), let[cel[ela,1]]), "FALSE")
      mk2 = length(x2)
      cltp2 = cl + 2
      y2 = udata[udata$clusters == br,]$AA.JUNCTION
      z2 = y2[x2]
      lengdif = FALSE
      
      # Find wich new sub-cluster has more sequences and give it first cluster name
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
    }else{
      # Find sequences contains the cell's similarity group in cell's position 
      strings.to.find = unlist(sim[cel[ela,1]])
      x1 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), str_c(strings.to.find, collapse="|")), "TRUE")
      mk1 = length(x1)
      cltp1 = cl + 1
      y1 = udata[udata$clusters == br,]$AA.JUNCTION
      z1 = y1[x1]
      
      # The other sequences of the cluster
      x2 = str_which(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), str_c(strings.to.find, collapse="|")), "FALSE")
      mk2 = length(x2)
      cltp2 = cl + 2
      y2 = udata[udata$clusters == br,]$AA.JUNCTION
      z2 = y2[x2]
      lengdif = FALSE
      
      # Find wich new sub-cluster has more sequences and give it first cluster name
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
    nn = FALSE # Return nn in it's original value 
    leaf2 = FALSE
  }
  
  ############################################ Old Control #############################################
  br = br + 1 # Increase the branch by 1
  met = met + 2 # Increase the counter by 2
  list2 = list("listax" = listax,"ggdf" = ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = ela, "cel" = cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
  if( ((clep[br-1] < clep[br]) && (sum(clep == clep[br]) != (2^clep[br]))) == TRUE && (length(which(udata$clusters == br)) > 1) ){ # If the next branch is in the next level
    met = geomSeq(1,2,1,1000)[ep+1]
  }
  while (length(which(udata$clusters == br)) <= 1){ # While the number of sequences in the branch is less than 2, go to the next branch and change counter 
    # Condition wich ends the algorithm 
    if(is.na(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]))){
      clmax = cl
      brtemp = br
      progend = TRUE
      list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
      if(brtemp < clmax){
        # Find the statistics of the remaining clusters before you end
        for (i in (br+1):clmax) {
          brtemp = i
          list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
          list2 = Matrices(list2,leaf,let,sim,d,algo,algocol,backcol,backcolj6,logFile)
        }
      }
      en = toc(quiet = TRUE)
      cat(paste0("Choice","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
      return(list2)
    }
    
    # Run Matrices as leaf
    leaf = TRUE
    list2 = list("listax" = listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = ela, "cel" = cel,"endper" = list2$endper, "last" = list2$last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
    list2 = Matrices(list2,leaf,let,sim,d,algo,algocol,backcol,backcolj6,logFile)
    
    # Condition wich ends the algorithm 
    if(is.na(((clep[br] < clep[br+1]) && (sum(clep == clep[br]) != (2^clep[br]))))){
      clmax = cl
      brtemp = br
      progend = TRUE
      list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
      if(brtemp < clmax){
        # Find the statistics of the remaining clusters before you end
        for (i in (br+1):clmax) {
          brtemp = i
          list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last,"progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
          list2 = Matrices(list2,leaf,let,sim,d,algo,algocol,backcol,backcolj6,logFile)
        }
      }
      en = toc(quiet = TRUE)
      cat(paste0("Choice","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
      return(list2)
    }
    
    # If the next branch is in the next level
    if( ((clep[br] < clep[br+1]) && (sum(clep == clep[br]) != (2^clep[br]))) == TRUE ){
      met = 0
      ep = ep + 1
      br = br + 1
    }else{
      br = br + 1
      met = met +2 
    }
    list2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = list2$last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
  }
  
  # When the counter reaches the end value (geometric sequence) we increase the level counter
  if( met == geomSeq(1,2,1,1000)[ep+1]){
    met = 0
    ep = ep + 1
  }
  
  # The return list
  leaf = FALSE
  result2 = list("listax" = list2$listax,"ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = udata,"permat"= list2$permat, "persim" = list2$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = ela, "cel" = cel,"endper" = list2$endper, "last" = list2$last, "progend" = progend, "leaf" = leaf,"leaf2" = leaf2)
  en = toc(quiet = TRUE)
  cat(paste0("Choice","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  return(result2)
}


######################################## Geometric Sequence ########################################
# A function that generates a geometric sequence
geomSeq <- function(start,ratio,begin,end){
  begin=begin-1
  end=end-1
  start*ratio**(begin:end)
}


######################################## Visualization 1 - Static Tree ########################################
# Creating a Static Tree visualization
Den <- function(lev,df,lastlist,levcut,Clus,flagtic,algo,threshold1,threshold2,enthr,logFile){
  if(flagtic == TRUE) tic()
  levcut = levcut +1
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
      }
    }
  }
  
  xN <- as.Node(x)
  if(flagtic == TRUE){
    en = toc(quiet = TRUE)
    cat(paste0(sprintf("Den -- Level: %d ", lev),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  }
  plot(xN)
}


######################################## Visualization 2a - SatLev ########################################
# Copmute the sequnces, Identity, Similarity of level and leaf (until this level) clusters and make barplots
SatLev <- function(lev,lastlist,Clus,df,flagtic,logFile){
  if(flagtic == TRUE) tic()
  if( is.element(lev,lastlist$clep) == FALSE){
    NULL
  }else{
    xN <- as.Node(df)
    # Level clusters
    t1 = which(lastlist$clep == lev)
    
    # leaves until level
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(lastlist$clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    
    # Visualization with barplot
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

######################################## Visualization 2b - LogoLev ########################################
# Create custom colour scheme
cs1 = make_col_scheme(chars=c("F","W","A","I","L","V","M","C","P","G","Y","T","S","H","K","R","E","D","Q","N"),
                      cols=c("#1E90FF", "#BA55D3", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#C6E2FF", "#C6E2FF", "#FFD700", "#00EE00", "#C1FFC1", "#54FF9F", "#54FF9F", "#FF0000", "#FF0000", "#FF0000", "#FFD700", "#FFD700", "#ED9121", "#ED9121"))

# A function to plot logos of level and leaf (until this level) clusters
LogoLev <- function(lev,lastlist,Clus,df,flagtic,seqen,ideen,simen,seqthr,idethr,simthr,logFile){
  if(flagtic == TRUE) tic()
  if( is.element(lev,lastlist$clep) == FALSE){
    if(flagtic == TRUE){
      en = toc(quiet = TRUE)
      cat(paste0(sprintf("LogoLev -- Level: %d ", lev),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    }
    NULL
  }else{
    xN <- as.Node(df)
    # Level clusters
    t1 = which(lastlist$clep == lev)
    
    # leaves until level
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(lastlist$clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    
    # Actions with Sequence Threshold Enabled
    if(seqen == TRUE){
      if(length(which(Clus[t1+1,]$seqnum < seqthr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$seqnum < seqthr)]
      }
    }
    
    # Actions with Identity Threshold Enabled
    if(ideen == TRUE){
      if(length(which(Clus[t1+1,]$Identity < idethr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$Identity < idethr)]
      }
    }
    
    # Actions with Similarity Threshold Enabled
    if (simen == TRUE){
      if(length(which(Clus[t1+1,]$Similarity <  simthr)) != 0){
        t1 = t1[-which(Clus[t1+1,]$Similarity <  simthr)]
      }
    }

    listff = list()
    nc = 3 # 3 logos in every row
    if(length(t1) > 0){
      for(i in 1:length(t1)){
        # leaves
        if(t1[i]<= max(xm2)){
          listff$temp = na.omit(df[df[names(df) == sprintf('level.%d', lastlist$clep[t1[i]])] == t1[i],]$AA.JUNCTION)
          names(listff)[length(listff)] = sprintf('Cluster.%d - seqnum:%d - Id:%g - Sim:%g - leaf', t1[i],Clus$seqnum[t1[i]+1],Clus$Identity[t1[i]+1],Clus$Similarity[t1[i]+1]) 
        }else{
          #level clusters
          x1 = as.data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", lev))] == t1[i], ]$AA.JUNCTION))
          if (nrow(x1) >0){
            names(x1)[1]= "AA.JUNCTION"
            x1 = as.character(x1$AA.JUNCTION)
            listff$temp = x1
            names(listff)[length(listff)] = sprintf('Cluster.%d - seqnum:%d - Id:%g - Sim:%g', t1[i],Clus$seqnum[t1[i]+1],Clus$Identity[t1[i]+1],Clus$Similarity[t1[i]+1]) 
          } 
        }
      }
      if(flagtic == TRUE){
        en = toc(quiet = TRUE)
        cat(paste0(sprintf("LogoLev -- Level: %d ", lev),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
      }
      # Plot logos
      ggseqlogo(listff, ncol=nc, method = "prob",col_scheme=cs1) 
    }else{
      if(flagtic == TRUE){
        en = toc(quiet = TRUE)
        cat(paste0(sprintf("LogoLev -- Level: %d ", lev),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
      }
      NULL
    }
  }
}


######################################## Visualization 2c - LogoCl ########################################
# A function to plot the logo of a cluster
LogoCl <- function(cl,lastlist,df,flagtic,logFile){
  if(flagtic == TRUE) tic()
  if(flagtic == TRUE){
    en = toc(quiet = TRUE)
    cat(paste0(sprintf("LogoCl -- Cluster: %d ", cl),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  }
  if(is.na(lastlist$clep[cl])){
   # if(flagtic == TRUE){
  #    en = toc(quiet = TRUE)
  #    cat(paste0(sprintf("LogoCl -- Cluster: %d ", cl),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  #  }
   NULL
  }else{
    if(flagtic == TRUE){
      en = toc(quiet = TRUE)
      cat(paste0(sprintf("LogoCl -- Cluster: %d ", cl),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    }
    xN <- as.Node(df)
    listff = list()
    listff$temp = na.omit(df[df[names(df) == sprintf('level.%d', lastlist$clep[cl])] == cl,]$AA.JUNCTION)
    xm = as.numeric(as.data.frame(xN$leaves))
    
    if(is.element(cl,xm)){
      names(listff)= sprintf('Cluster: %s    num:%d - leaf', cl ,lastlist$dfsum$len[cl+1])  
    }else{
      names(listff)= sprintf('Cluster: %s    num:%d', cl ,lastlist$dfsum$len[cl+1]) 
    }
    # Plot logo
    ggseqlogo(listff, method = "prob", col_scheme=cs1)
  }
}

######################################## Visualization 3a - BarLev ########################################
# A function to plot barplots of level and leaf (until this level) clusters
BarLev <- function(lev,lastlist,df,perlist,persimlist,Clus,let,sim,cho,flagtic,logFile){
  if(flagtic == TRUE) tic()
  if(cho == "Identity"){
    if( is.element(lev,lastlist$clep) == FALSE){
      plot(1, type="n", axes=F, xlab="", ylab="")
    }else{
      xN <- as.Node(df)
      # Level clusters
      t2 = which(lastlist$clep == lev)
      
      # leaves until level
      xm = as.numeric(as.data.frame(xN$leaves))
      orio = min(which(lastlist$clep == lev))
      xm2 = sort(xm[xm < orio])
      t2 = sort(append(xm2,t2,after = length(xm2)))
      
      # 3 barplots in every row
      if(length(t2) %% 3 == 0){
        par(mfrow = c(length(t2)%/%3,3))
      }else{
        par(mfrow = c(length(t2)%/%3 + 1,3))
      }
      for(i in 1:length(t2)){
        ar = str_which(names(perlist),as.character(t2[i]))[1]
        par(xpd=TRUE)
        output <- matrix(unlist(perlist[ar]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
        if(i<= length(xm2)){ # leaves
          barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d - seqnum:%d - leaf', t2[i],Clus$seqnum[t2[i]+1])) 
          
        }else{ # level clusters
          barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d - seqnum:%d', t2[i],Clus$seqnum[t2[i]+1])) 
        }
        legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=let,cex = 0.6)
      }
    }
  }else{
    if( is.element(lev,lastlist$clep) == FALSE){
      plot(1, type="n", axes=F, xlab="", ylab="")
    }else{
      xN <- as.Node(df)
      # Level clusters
      t2 = which(lastlist$clep == lev)
      
      # leaves until level
      xm = as.numeric(as.data.frame(xN$leaves))
      orio = min(which(lastlist$clep == lev))
      xm2 = sort(xm[xm < orio])
      t2 = sort(append(xm2,t2,after = length(xm2)))
      
      # 3 barplots in every row
      if(length(t2) %% 3 == 0){
        par(mfrow = c(length(t2)%/%3,3))
      }else{
        par(mfrow = c(length(t2)%/%3 + 1,3))
      }
      for(i in 1:length(t2)){
        ar = str_which(names(persimlist),as.character(t2[i]))[1]
        par(xpd=TRUE)
        output <- matrix(unlist(persimlist[ar]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
        if(i<= length(xm2)){ # leaves
          barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d - seqnum:%d - leaf', t2[i],Clus$seqnum[t2[i]+1])) 
          
        }else{ # level clusters
          barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d - seqnum:%d', t2[i],Clus$seqnum[t2[i]+1])) 
        }
        legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=names(sim),cex = 0.6)
      }
    }
  }
  if(flagtic == TRUE){
    en = toc(quiet = TRUE)
    cat(paste0(sprintf("BarLev -- Level: %d ", lev),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  }
}


######################################## Visualization 3b - BarCl ########################################
# A function to plot the barplot of a cluster
BarCl <- function(cl,perlist,persimlist,Clus,let,sim,cho,lastlist,flagtic,logFile){
  if(flagtic == TRUE) tic()
  if(cho == "Identity"){
    if(is.na(lastlist$clep[cl])){
      plot(1, type="n", axes=F, xlab="", ylab="")
    }else{
      par(mar=c(3,3,4,4),xpd=TRUE)
      output <- matrix(unlist(perlist[cl+1]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
      barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d',cl)) 
      legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=let,cex = 0.6)
    }
  }else{
    if(is.na(lastlist$clep[cl])){
      plot(1, type="n", axes=F, xlab="", ylab="")
    }else{
      par(mar=c(3,3,4,4),xpd=TRUE)
      output <- matrix(unlist(persimlist[cl+1]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
      barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d',cl)) 
      legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=names(sim),cex = 0.6) 
    }
  }
  if(flagtic == TRUE){
    en = toc(quiet = TRUE)
    cat(paste0(sprintf("BarCl -- Cluster :%d ", cl),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  }
}


######################################## Visualization 4a - AminoLev - Unused ########################################
# A function for showing  Sequences and Id's for specific level
AminoLev <- function(level,lastlist,df,Clus,flagtic,logFile){
  if(flagtic == TRUE) tic()
  if( is.element(level,lastlist$clep) == FALSE){
    if(flagtic == TRUE){
      en = toc(quiet = TRUE)
      cat(paste0(sprintf("AminoLev -- Level: %d ", level),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    }
    x4 = table(NA)
  }else{
    # Sequences of Level
    sum(Clus[Clus$level == level,]$seqnum)
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
    if(flagtic == TRUE){
      en = toc(quiet = TRUE)
      cat(paste0(sprintf("AminoLev -- Level: %d ", level),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    }
    x4 # Print in console 
  }
}


######################################## Visualization 4b - AminoCl ########################################
# A function for showing  Sequences and Id's for specific cluster
AminoCl <- function(clust,lastlist,df,flagtic,logFile){
  #if(flagtic == TRUE) tic("AminoCl")
  #if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  if(clust == 0){
    x3 = data.frame("Sequence.ID" = df$Sequence.ID, "AA.JUNCTION" = df$AA.JUNCTION, "V.GENE.and.allele" = df$V.GENE.and.allele)
    x3
  }else{
    if(is.na(lastlist$clep[clust]) == FALSE){
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


######################################## Visualization 5a - Opt (Absolute Consensus for letters) ########################################
# A function which visualize the sequences of a data frame using absolute common letters (i.e. "A _ _ _ _ K R _ _ _ Q _ Y Y Y _ _ _ T _")
Opt <- function(df,flagtic,logFile){
  if(is.list(df)){
    #if(flagtic == TRUE) tic()
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
      xar[i] = "_"
      if (f[1] == nrow(df)){
        xar[i] = names(f[1])
      }
      
    }
    #if(flagtic == TRUE){
    #  en = toc(quiet = TRUE)
    #  cat(paste0("Opt","\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    #}
    # String Version
    #paste(xar,collapse = ' ')
    # Matrix Version
    xar
  }
}


######################################## Visualization 5b - OptNew (Percentage Consensus for letters) ########################################
# A function which visualize the sequences and percentages of a data frame using percentage common letters (i.e. "A _ _ _ _ K R _ _ _ Q _ Y Y Y _ _ _ T _")
OptNew <- function(df,logFile){
  if(is.list(df)){
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    xarpos <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
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
    rbind(xar,xarpos)
  }
}


######################################## Visualization 5c - Opt2 (Absolute Consensus for similarity groups) ########################################
# A function which visualize the sequences of a data frame using absolute similarity groups (i.e. "A _ _ _ _ K Am _ Ba _ Ac _ Y Y Y _ _ _ T _")
Opt2 <- function(df,sim,alt,altsim,flagtic,logFile){
  if(is.list(df)){
    #if(flagtic == TRUE) tic()
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
      xar[i] = "_"
      if (f[1] == nrow(df)){
        xar[i] = names(f[1])
      }else{
        y = TRUE
        d = str_which(names(sim),gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),names(f[1])))])))
        for(j in 2:length(f)){
          y = y && is.element(names(f[j]),unlist(sim[d]))
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
    # String Version
    #paste(xar,collapse = ' ')
    # Matrix Version
    xar
  }
}

######################################## Visualization 5d - Opt2New (Percentage Consensus for similarity groups) ########################################
# A function which visualize the sequences and percentages of a data frame using absolute similarity groups (i.e. "A _ _ _ _ K Am _ Ba _ Ac _ Y Y Y _ _ _ T _")
Opt2New <- function(df,sim,logFile){
  if(is.list(df)){
    xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    xarpos <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
    for (i in 1:str_length(df$AA.JUNCTION[1])) {
      f <- table(str_sub(df$AA.JUNCTION,i,i))
      for(j in 1:length(f)){
        if( length(which(str_detect(names(f),names(sim[str_which(names(sim),gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),names(f[j])))])))])))) == 0 ){
          names(f)[j] = gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),names(f[j])))]))
        }else if(length(which(str_detect(names(f),names(sim[str_which(names(sim),gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),names(f[j])))])))])))) != j){
          f[which(str_detect(names(f),names(sim[str_which(names(sim),gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),names(f[j])))])))])))] = f[which(str_detect(names(f),names(sim[str_which(names(sim),gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),names(f[j])))])))])))] + f[j]
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
    #if(flagtic == TRUE){
    #  en = toc(quiet = TRUE)
    #  cat(paste0("Opt2New","\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    #}
    rbind(xar,xarpos)
  }
}


######################################## Visualization 6 - Id (Plots of Identity or Similarity) ########################################
# A function plot a diagramm with Identity or Similarity of clusters in every level
Id <- function(ff,cho,flagtic,logFile){
  if(flagtic == TRUE) tic()
  if(cho == "Identity"){
    par(xpd=TRUE)
    pp <- data.frame(x = ff$level,y = ff$sumper)
    if(flagtic == TRUE){
      en = toc(quiet = TRUE)
      cat(paste0(sprintf("Id -- Type: $s ", cho),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    }
    ggplot(na.omit(pp),aes(x = x,y = y)) + stat_sum()
  }else{
    par(xpd=TRUE)
    pp <- data.frame(x = ff$level,y = ff$sumper2)
    if(flagtic == TRUE){
      en = toc(quiet = TRUE)
      cat(paste0(sprintf("Id -- Type: $s ", cho),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    }
    ggplot(na.omit(pp),aes(x = x,y = y)) + stat_sum()
  }
}


######################################## Visualization 7 - cc (Collapsible Tree) ########################################
# A function plots a collapseble tree with id, sequences, identity and similarity of every cluster - node 
cc <- function(df,Clus,flagtic,logFile){
  if(flagtic == TRUE) tic()
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
  if(flagtic == TRUE){
    en = toc(quiet = TRUE)
    cat(paste0("Collapsible Tree","\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  }
  collapsibleTree(
    trid1,
    hierarchy = colnames(df[str_which(names(df), "level.")]),
    fill = c("jj",na.omit(ii)),
    width = 1820,
    height = 775,
    collapsed = FALSE
  )
}


######################################## Visualization 8a - idenlev (Identity or Similarity Values) ########################################
# A function show the value of identity or similarity for level clusters and leaves until this level
idenlev <- function(lev,Clus,df,cho,flagtic,logFile){
  if(flagtic == TRUE) tic()
  if(cho == "Identity"){
    if(flagtic == TRUE){
      en = toc(quiet = TRUE)
      cat(paste0(sprintf("idenlev -- Level: %d ", lev),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    }
    xN <- as.Node(df)
    # Level clusters
    t1 = which(lastlist$clep == lev)
    
    # leaves until level
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(lastlist$clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    sum(na.omit(Clus[t1+1,]$Identity)) /  length(na.omit(Clus[t1+1,]$Identity))
  }else{
    if(flagtic == TRUE){
      en = toc(quiet = TRUE)
      cat(paste0(sprintf("idenlev -- Level: %d ", lev),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    }
    # Level clusters
    t1 = which(lastlist$clep == lev)
    
    # leaves until level
    xm = as.numeric(as.data.frame(xN$leaves))
    orio = min(which(lastlist$clep == lev))
    xm2 = sort(xm[xm < orio])
    t1 = sort(append(xm2,t1,after = length(xm2)))
    sum(na.omit(Clus[t1+1,]$Similarity)) /  length(na.omit(Clus[t1+1,]$Similarity))
  }
}


######################################## Visualization 8b - idencl (Identity or Similarity Values) ########################################
# A function show the value of identity or similarity for a cluster
idencl <- function(cl,Clus,cho,flagtic,logFile){
  if(flagtic == TRUE) tic()
  if(is.element(cl,Clus$ClusterId)){
    if(cho == "Identity"){
      if(flagtic == TRUE){
        en = toc(quiet = TRUE)
        cat(paste0(sprintf("idencl -- Cluster :%d ", cl),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
      }
      na.omit(Clus[Clus$ClusterId == cl,])$Identity
    }else{
      if(flagtic == TRUE){
        en = toc(quiet = TRUE)
        cat(paste0(sprintf("idencl -- Cluster :%d ", cl),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
      }
      na.omit(Clus[Clus$ClusterId == cl,])$Similarity 
    }
    if(flagtic == TRUE){
      en = toc(quiet = TRUE)
      cat(paste0(sprintf("idencl -- Cluster :%d ", cl),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
    }
  }
}


######################################## Visualization 8c - EmPin (Identity or Similarity Values) ########################################
# A function create a table with average identity, identity standard deviation, average similarity and similarity standard deviation for every level 
EmPin <- function(lastlist,Clus,df,dfsd,flagtic,logFile){
  if(flagtic == TRUE) tic()
  xN <- as.Node(df)
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
  if(flagtic == TRUE){
    en = toc(quiet = TRUE)
    cat(paste0("EmPin","\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  }
  dfsd
}


######################################## Visualization 9 - Netw (Network) ########################################
# A function creates a network for a given level 
Netw <- function(lev,thr,thrt,netyp,df,lastlist,Clus,sim,altsim,ts,shth,net_sil,flagtic,algo,threshold1,threshold2,enthr,logFile,levcut){
  if(flagtic == TRUE) tic()
  
  # Compute tree
  levcut = levcut + 1
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
  
  # If threshold is enabled 
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
      }
    }
    xN <- as.Node(x)
    
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
  
  # Clusters of given level
  bo = which(lastlist$clep == lev)
  ffg = as.data.frame(matrix(0,nrow = max(bo)+1,ncol = (max(bo)+1))) # Matrix with the distance between 2 Nodes
  ffg2 = as.data.frame(matrix(0,nrow = max(bo)+1,ncol = (max(bo)+1))) # Matrix with letter similarities between 2 Nodes 
  ffg2sim = as.data.frame(matrix(0,nrow = max(bo)+1,ncol = (max(bo)+1))) # Matrix with group similarities between 2 Nodes 
  ffg3 = as.data.frame(matrix(NA,nrow = max(bo)+1,ncol = (max(bo)+1))) # the final distance combining ffg and ffg2 and fgg2sim
  
  if(enthr == TRUE){
    bfbf = unique(as.numeric(unlist(sss)))
    # Computation of Matrices
    s = sapply(1:length(bfbf), function(i){
      hhh = str_which(sss,sprintf("%d", bfbf[i]))
      tempN1 = unlist(sss[hhh[1]])
      # Distance of Node 1 from root of the tree
      tempN1 = tempN1[1:str_which(tempN1,sprintf("%d", bfbf[i]))]
      # Consensus for letter of Node 1 
      temp2N1 = paste(Opt(AminoCl(bfbf[i],lastlist,df,flagtic,logFile),flagtic,logFile),collapse = ' ')
      # Consensus for groups Similarity of Node 1 
      temp3N1 = paste(Opt2(AminoCl(bfbf[i],lastlist,df,flagtic,logFile),sim,TRUE,altsim,flagtic,logFile),collapse = ' ')
      sapply(i:length(bfbf),function(j) {
        hhh2 = str_which(sss,sprintf("%d", bfbf[j]))
        tempN2 = unlist(sss[hhh2[1]])
        # Distance of Node 2 from root of the tree
        tempN2 = tempN2[1:str_which(tempN2,sprintf("%d", bfbf[j]))]
        tem = which(tempN1 %in% tempN2)
        d1 = length(tempN1) - tem[length(tem)]
        d2 = length(tempN2) - tem[length(tem)]
        temp2N2 = paste(Opt(AminoCl( bfbf[j],lastlist,df,flagtic,logFile),flagtic,logFile),collapse = ' ')
        # Consensus for letter of Node 2
        temp2N2 = str_replace_all(temp2N2,"_"," ")
        temp3N2 = paste(Opt2(AminoCl( bfbf[j],lastlist,df,flagtic,logFile),sim,TRUE,altsim,flagtic,logFile),collapse = ' ')
        # Consensus for groups Similarity of Node 2
        temp3N2 = str_replace_all(temp3N2,"_"," ")
        # Set values in matrices cell
        set(ffg,  bfbf[i]+1,  bfbf[j]+1, d1 + d2) 
        set(ffg2, (bfbf[i]+1),  bfbf[j]+1, stringdist(temp2N1,temp2N2,method = "lv" )) 
        set(ffg2sim, (bfbf[i]+1),  bfbf[j]+1, stringdist(temp3N1,temp3N2,method = "lv" )) 
      })
    })
  }else{
    # Computation of Matrices
    s = sapply(1:(max(bo)+1), function(i){
      tempN1 = FindNode(xN,(sprintf("%d", i-1)))
      # Distance of Node 1 from root of the tree
      tempN1 = tempN1$path
      # Consensus for letter of Node 1 
      temp2N1 = paste(Opt(AminoCl(i-1,lastlist,df,flagtic,logFile),flagtic,logFile),collapse = ' ')
      # Consensus for groups Similarity of Node 1 
      temp3N1 =  paste(Opt2(AminoCl(i-1,lastlist,df,flagtic,logFile),sim,TRUE,altsim,flagtic,logFile),collapse = ' ')
      sapply(i:(max(bo)+1),function(j) {
        tempN2 = FindNode(xN,(sprintf("%d", j-1)))
        # Distance of Node 2 from root of the tree
        tempN2 = tempN2$path
        tem = which(tempN1 %in% tempN2)
        d1 = length(tempN1) - tem[length(tem)]
        d2 = length(tempN2) - tem[length(tem)]
        temp2N2 = paste(Opt(AminoCl(j-1,lastlist,df,flagtic,logFile),flagtic,logFile),collapse = ' ')
        # Consensus for letter of Node 2
        temp2N2 = str_replace_all(temp2N2,"_"," ")
        temp3N2 = paste(Opt2(AminoCl(j-1,lastlist,df,flagtic,logFile),sim,TRUE,altsim,flagtic,logFile),collapse = ' ')
        # Consensus for groups Similarity of Node 2
        temp3N2 = str_replace_all(temp3N2,"_"," ")
        # Set values in matrices cell
        set(ffg, i, j, as.integer(d1 + d2)) 
        set(ffg2, i, j, stringdist(temp2N1,temp2N2,method = "lv" )) 
        set(ffg2sim, i, j, stringdist(temp3N1,temp3N2,method = "lv" )) 
      })
    })
  }
  
  #find max
  ma = max(ffg[is.na(ffg) == FALSE])
  ffg[ffg == 0] = ma
  ffg2[ffg2 == 0] = str_length(lastlist$udata$AA.JUNCTION[1])
  ffg2sim[ffg2sim == 0] = str_length(lastlist$udata$AA.JUNCTION[1])
  ffg[lower.tri(ffg)] = NA
  ffg2[lower.tri(ffg2)] = NA
  ffg2sim[lower.tri(ffg2sim)] = NA
  ffg = ffg / ma * 100
  ffg2 = str_length(lastlist$udata$AA.JUNCTION[1]) - ffg2
  ffg2 = ffg2 / str_length(lastlist$udata$AA.JUNCTION[1]) * 100
  ffg2sim = str_length(lastlist$udata$AA.JUNCTION[1]) - ffg2sim
  ffg2sim = ffg2sim / str_length(lastlist$udata$AA.JUNCTION[1]) * 100
  ffg3 = (2/str_length(lastlist$udata$AA.JUNCTION[1])) * ffg2 + (2/str_length(lastlist$udata$AA.JUNCTION[1])) * ffg2sim + (2 / ma) * (ma - ffg)
  
  # Not show self connenctions 
  diag(ffg) = 0
  diag(ffg2) = 0
  diag(ffg2sim) = 0
  diag(ffg3) = 0
  ffg3 = (ffg2 + ffg2sim + ffg) / 3
  
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
  tempor[thrtyp > thr] = 0
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
  
  #if silhouette enabled
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
  
  # Create matrix for ploting the network
  jjj = matrix(1,nrow = max(bo)+1,ncol = (max(bo)+1))
  net0 <- graph_from_adjacency_matrix(jjj,mode = "upper")
  net1 <- graph_from_adjacency_matrix(jjj,mode = "upper")
  
  bb = vector(length = (max(bo)+1))
  bb[1]=0
  bb[2:length(bb)] = lastlist$clep[1:max(bo)]
  #XRWMA
  # Generate colors based on media type:
  colrs <- c("#1E90FF", "#BA55D3", "#0000FF", "#557fd2", "#54d17e", "#8aad62", "#C6E2FF", "#e5e234", "#FFD700", "#00EE00", "#C1FFC1", "#ea8509", "#54FF9F", "#FF0000", "#ed3b1c", "#ed1c7a", "#0c0c0c", "#b8d8af", "#ED9121","#45f713")
  colrs = colrs[1:(max(na.omit(bb))+1)]
  V(net0)$color <- colrs[bb+1]
  V(net1)$color <- colrs[bb+1]
  
  # Compute node degrees (#links) and use that to set node size:
  deg <- (Clus$seqnum[1:length(bb)] / Clus$seqnum[1]) * 20
  V(net0)$size <- deg
  V(net1)$size <- deg
  # The labels are currently node IDs.
  # Setting them to NA will render no labels:
  V(net0)$label <- Clus$ClusterId[1:length(bb)]
  V(net1)$label <- Clus$ClusterId[1:length(bb)]
  
  # Set edge width based on weight:
  hhh = na.omit(as.vector(t(ffg3)))
  hhh1 = na.omit(as.vector(t(tempor3)))
  E(net0)$width <- hhh * 0.06 
  E(net1)$width <- hhh1 
  
  #change arrow size and edge color:
  E(net0)$arrow.size <- .2
  E(net0)$edge.color <- "gray80"
  E(net1)$arrow.size <- .2
  E(net1)$edge.color <- "gray80"
  pal1 <- rainbow(6, alpha=1) 
  
  net0.copy <- igraph::delete.edges(net0, which(E(net0)$width == 0))
  net0.copy <- igraph::delete.vertices(net0.copy, is.na(V(net0)$size))
  net1.copy <- igraph::delete.edges(net1, which(E(net1)$width == 0))
  
  if(flagtic == TRUE){
    en = toc(quiet = TRUE)
    cat(paste0(sprintf("Network -- level: %d, Network type: %s, Threshold type: %s, Threshold value: %d, Using silhouette: %d ", lev,netyp,thrt,thr,net_sil),"\t","-","\t","-","\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  }
  listnet = list("net0.copy" = net0.copy, "net1.copy" = net1.copy, "bb" = bb, "colrs" = colrs, "ma" = ma, "ffg" = ffg, "ffg2" = ffg2, "ffg2sim" = ffg2sim, "ffg3" = ffg3 ,"tempor" = tempor, "tempor2" = tempor2, "tempor2sim" = tempor2sim, "tempor3" = tempor3,"thrtyp" = thrtyp, "thrtyp2" = thrtyp2)
  return(listnet)
}
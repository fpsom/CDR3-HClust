# Function Matrices compute the apsolute matrix, the matrix with percentages, the absolute similarity matrix and the percentage similarity matrix
Matrices <- function(list1,leaf,let,sim,d,algo,algocol){
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
  print(udata)
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
  result1=list("ggdf" = list1$ggdf, "dfsum" = list1$dfsum,"list" = listxx,"listn" = listyy, "udata" = udata,"permat"= permat, "persim" = persim, "br" = br, "cl" = cl, "met" = list1$met, "ep" = list1$ep, "clep" = list1$clep, "nn" = list1$nn, "sumper" =list1$sumper,"sumper2" = list1$sumper2, "ela" = list1$ela, "cel" = list1$cel, "endper" = list1$endper, "last" = list1$last)
  if ( leaf == TRUE){
    Finish(result1,TRUE,let,sim,d,algo,algocol)
  }else{
    Finish(result1,FALSE,let,sim,d,algo,algocol)
  }
}


Finish <- function(list2,leaf,let,sim,d,algo,algocol){
  udata = list2$udata
  permat = list2$permat
  persim = list2$persim
  sumper = list2$sumper # The total percentage of permat 
  sumper2 = list2$sumper2 # The total percentage of persim
  endper = list2$endper
  dfsum = list2$dfsum
  br = list2$br
  cl = list2$cl
  nn = list2$nn
  last = list2$last
  t1 =which(permat[-length(permat),] == 100,arr.ind = TRUE)
  sumper = (length(as.numeric(t1[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
  t2 =which(persim[-length(persim),] == 100,arr.ind = TRUE)
  sumper2 = (length(as.numeric(t2[,2]))* 100) / str_length(udata$AA.JUNCTION[1])
  if(algo == "Identity"){
    if (sumper > endper){
      nn = TRUE  # When nn = TRUE the percentage of sumper < endper%
      if (sumper2 > endper){
        last = last+1
        if (last == 1 && leaf == FALSE){
          clmax = cl
          brtemp = br
          while (brtemp < clmax) {
            brtemp = brtemp +1 
            print(brtemp)
            print(clmax)
            print(last)
            list2 = list("ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = last)
            list2 = Matrices(list2,TRUE,let,sim,d,algo,algocol)
          }
        }
        return(list2)  # End of the programm, returning the final list
      }
    }
  }else{
    if (sumper2 > endper){
      last = last+1
      if (last == 1 && leaf == FALSE){
        clmax = cl
        brtemp = br
        while (brtemp < clmax) {
          brtemp = brtemp +1 
          print(brtemp)
          print(clmax)
          print(last)
          list2 = list("ggdf" = list2$ggdf, "dfsum" = list2$dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = brtemp, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = list2$nn, "sumper" = list2$sumper, "sumper2" = list2$sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = last)
          list2 = Matrices(list2,TRUE,let,sim,d,algo,algocol)
        }
      }
      return(list2)  # End of the programm, returning the final list
    }
  }
  vv = length(udata[udata$clusters == br,]$AA.JUNCTION)
  dfsum[nrow(dfsum) + 1,] = c(sumper,sumper2,br,vv)
  result2 = list("ggdf" = list2$ggdf, "dfsum" = dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = list2$br, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = last)
  if ( leaf == TRUE){
    nn = FALSE
    result2 = list("ggdf" = list2$ggdf, "dfsum" = dfsum,"list" = list2$list, "listn" = list2$listn, "udata" = list2$udata,"permat"= list2$permat, "persim" = list2$persim, "br" = list2$br, "cl" = list2$cl, "met" = list2$met, "ep"= list2$ep, "clep" = list2$clep,"nn" = nn, "sumper" = sumper, "sumper2" = sumper2, "ela" = list2$ela, "cel" = list2$cel,"endper" = list2$endper, "last" = last)
    return(result2)
  }else{
    Choice(result2,let,sim,d,algo,algocol)
  }
}


# Function Choice choose which matrix cell will be used for the division of the data 
Choice <- function(list3,let,sim,d,algo,algocol){
  nn = list3$nn
  br = list3$br
  cl = list3$cl
  if(nn == TRUE || (algo == "Similarity") ){
    permat = list3$persim # If sumper < endper% we want to check only the persim matrix
  }else{
    permat =list3$permat  # Else permat and if it is necessary the persim matrix
    persim = list3$persim
  }
  cel = which(permat[,(algocol + 1):ncol(permat)] == max(permat[,(algocol + 1):ncol(permat)]), arr.ind = TRUE)
  ela = 1
  poss = max(permat[,(algocol + 1):ncol(permat)]) 
  if (max(permat[,(algocol + 1):ncol(permat)]) == 100){ # We exclude the 100 % from the max values
    cel = which(permat[,(algocol + 1):ncol(permat)] == max(permat[,(algocol + 1):ncol(permat)][permat[,(algocol + 1):ncol(permat)]!=max(permat[,(algocol + 1):ncol(permat)])]), arr.ind = TRUE) # The desired cell
    poss = max(permat[,(algocol + 1):ncol(permat)][permat[,(algocol + 1):ncol(permat)]!=max(permat[,(algocol + 1):ncol(permat)])])
  }
  
  if(poss == 0){
    clmax = cl
    brtemp = br
    while (brtemp < clmax) {
      brtemp = brtemp +1 
      print(brtemp)
      list3 = list("ggdf" = list3$ggdf, "dfsum" = list3$dfsum,"list" = list3$list, "listn" = list3$listn, "udata" = list3$udata,"permat"= list3$permat, "persim" = list3$persim, "br" = brtemp, "cl" = list3$cl, "met" = list3$met, "ep"= list3$ep, "clep" = list3$clep,"nn" = list3$nn, "sumper" = list3$sumper, "sumper2" = list3$sumper2, "ela" = list3$ela, "cel" = list3$cel,"endper" = list3$endper, "last" = list3$last)
      list3 = Matrices(list3,TRUE,let,sim,d,algo,algocol)
    }
    return(list3)
  }
  
  if(algo == "Identity"){
    if ((length(cel)/2) > 1){
      eqcol = vector('numeric') # A vector for the cell's rows with the same entropy
      j = 1
      for(i in 1:(length(cel)/2)){
        if (permat[,(algocol + 1):ncol(permat)][nrow(permat[,(algocol + 1):ncol(permat)]),cel[i,2]] < permat[,(algocol + 1):ncol(permat)][nrow(permat[,(algocol + 1):ncol(permat)]),cel[ela,2]]) {  # If the percentage is the same (cel multidimensional) we keep the cel with the lowest entropy
          ela = i
          eqcol[] = 0
          eqcol[1] = i
          j = 2
        }else if(permat[,(algocol + 1):ncol(permat)][nrow(permat[,(algocol + 1):ncol(permat)]),cel[i,2]] == permat[nrow(permat[,(algocol + 1):ncol(permat)]),cel[ela,2]]){
          eqcol[j] = i
          j = j+1
        }
      }
      eqcol <- eqcol[-which(eqcol == 0)]
      if (length(eqcol)>1 && nn == FALSE){ # If the vector has 2 or more numbers means that we have columns with the same entropy and nn = FALSE in order not to double check the persim
        ela = eqcol[1]
        meg = persim[,(algocol + 1):ncol(persim)][str_which(sim,let[cel[ela,1]]),cel[ela,2]]
        for(i in 2:length(eqcol)){
          if(persim[,(algocol + 1):ncol(persim)][str_which(sim,let[cel[eqcol[i],1]]),cel[eqcol[i],2]] > meg){ # We find the maximum percentage of the similarity percentage marix
            meg = persim[,(algocol + 1):ncol(persim)][str_which(sim,let[cel[eqcol[i],1]]),cel[eqcol[i],2]]
            ela = eqcol[i]
          } else if (persim[,(algocol + 1):ncol(persim)][str_which(sim,let[cel[eqcol[i],1]]),cel[eqcol[i],2]] == meg){ # If the percentage is the same we use the entropy
            if(persim[,(algocol + 1):ncol(persim)][nrow(persim[,(algocol + 1):ncol(persim)]),cel[eqcol[i],2]] < persim[,(algocol + 1):ncol(persim)][nrow(persim[,(algocol + 1):ncol(persim)]),cel[ela,2]]){
              ela = i
            }
          }
        }
      }
    }
  }else{
    if ((length(cel)/2) > 1){
      ela = 1 
      for(i in 2:(length(cel)/2)){
        if (permat[,(algocol + 1):ncol(permat)][nrow(permat[,(algocol + 1):ncol(permat)]),cel[i,2]] < permat[,(algocol + 1):ncol(permat)][nrow(permat[,(algocol + 1):ncol(permat)]),cel[ela,2]]) {  # If the percentage is the same (cel multidimensional) we keep the cel with the lowest entropy
          ela = i
        }
      }
    }
  }
  nn = FALSE # Return nn in it's original value 
  result3 = list("ggdf" = list3$ggdf, "dfsum" = list3$dfsum,"list" = list3$list, "listn" = list3$listn, "udata" = list3$udata,"permat"= list3$permat, "persim" = list3$persim, "br" = list3$br, "cl" = list3$cl, "met" = list3$met, "ep"= list3$ep, "clep" = list3$clep,"nn" = nn, "sumper" = list3$sumper, "sumper2" = list3$sumper2, "ela" = ela, "cel" = cel,"endper" = list3$endper, "last" = list3$last)
  Divide(result3,let,sim,d,algo,algocol)
}


# Function Divide divide the data into 2 new clusters and updates the column with the right level 
Divide <- function(list4,let,sim,d,algo,algocol){
  ggdf = list4$ggdf
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
      ggdf[nrow(ggdf) + 1,] = c(cltp1,mk2)
      ggdf[nrow(ggdf) + 1,] = c(cltp2,mk1)
    }else{
      ggdf[nrow(ggdf) + 1,] = c(cltp1,mk1)
      ggdf[nrow(ggdf) + 1,] = c(cltp2,mk2)
    }
    for(i in 1:length(z1)){
      udata[str_which(udata$AA.JUNCTION,z1[i]),sprintf('level.%d', ep)] = cl+1
    }
    for(i in 1:length(z2)){
      udata[str_which(udata$AA.JUNCTION,z2[i]),sprintf('level.%d', ep)] = cl+2
    }
    if(lengdif == TRUE){
      udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), let[cel[ela,1]]), cl+2 ,cl+1)
    }else{
      udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), let[cel[ela,1]]), cl+1 ,cl+2)
    }
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
      ggdf[nrow(ggdf) + 1,] = c(cltp1,mk2)
      ggdf[nrow(ggdf) + 1,] = c(cltp2,mk1)
    }else{
      ggdf[nrow(ggdf) + 1,] = c(cltp1,mk1)
      ggdf[nrow(ggdf) + 1,] = c(cltp2,mk2)
    }
    for(i in 1:length(z1)){
      udata[str_which(udata$AA.JUNCTION,z1[i]),sprintf('level.%d', ep)] = cl+1
    }
    for(i in 1:length(z2)){
      udata[str_which(udata$AA.JUNCTION,z2[i]),sprintf('level.%d', ep)] = cl+2
    }
    if(lengdif == TRUE){
      udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), str_c(strings.to.find, collapse="|")), cl+2 ,cl+1) 
    }else{
      udata[udata$clusters == br,]$clusters <- ifelse(str_detect(str_sub(udata[udata$clusters == br,]$AA.JUNCTION,(cel[ela,2]+algocol),(cel[ela,2]+algocol)), str_c(strings.to.find, collapse="|")), cl+1 ,cl+2) 
    }
  }
  
  result4 = list("ggdf" = ggdf,"dfsum" = list4$dfsum,"list" = list4$list, "listn" = list4$listn, "udata" = udata,"permat"= list4$permat, "persim" = list4$persim, "br" = list4$br, "cl" = list4$cl, "met" = list4$met, "ep"= list4$ep, "clep" = list4$clep,"nn" = list4$nn, "sumper" = list4$sumper, "sumper2" = list4$sumper2, "ela" = list4$ela, "cel" = list4$cel,"endper" = list4$endper, "last" = list4$last)
  Control(result4,let,sim,d,algo,algocol)
}


# Function Control update the branch, the counter and the cluster numbers in order to begin a new valid sub-division
Control <- function(list5,let,sim,d,algo,algocol){
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
  list5 = list("ggdf" = list5$ggdf, "dfsum" = list5$dfsum,"list" = list5$list, "listn" = list5$listn, "udata" = list5$udata,"permat"= list5$permat, "persim" = list5$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list5$nn, "sumper" = list5$sumper, "sumper2" = list5$sumper2, "ela" = list5$ela, "cel" = list5$cel,"endper" = list5$endper, "last" = list5$last)
  if( ((clep[br-1] < clep[br]) && (sum(clep == clep[br]) != (2^clep[br]))) == TRUE && (length(which(udata$clusters == br)) > 2) ){ # If the next branch is in the next level
    met = geomSeq(1,2,1,50)[ep+1]
  }
  while (length(which(udata$clusters == br)) <= 2){ # While the number of sequences in the branch is less than 2, go to the next branch and change counter 
    if(is.na(str_length(udata[udata$clusters == br,]$AA.JUNCTION[1]))){
      return(list5)
    }
    list5 = Matrices(list5, TRUE,let,sim,d,algo,algocol)
    if(is.na(((clep[br] < clep[br+1]) && (sum(clep == clep[br]) != (2^clep[br]))))){
      return(list5)
    }
    if( ((clep[br] < clep[br+1]) && (sum(clep == clep[br]) != (2^clep[br]))) == TRUE ){ # If the next branch is in the next level
      met = 0
      ep = ep + 1
      br = br + 1
    }else{
      br = br + 1
      met = met +2 
    }
    list5 = list("ggdf" = list5$ggdf, "dfsum" = list5$dfsum,"list" = list5$list, "listn" = list5$listn, "udata" = list5$udata,"permat"= list5$permat, "persim" = list5$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list5$nn, "sumper" = list5$sumper, "sumper2" = list5$sumper2, "ela" = list5$ela, "cel" = list5$cel,"endper" = list5$endper, "last" = list5$last)
  }
  if( met == geomSeq(1,2,1,50)[ep+1]){ # When the counter reaches the end value (geometric sequence) we increase the level counter
    met = 0
    ep = ep + 1
  }
  result5 = list("ggdf" = list5$ggdf, "dfsum" = list5$dfsum,"list" = list5$list, "listn" = list5$listn, "udata" = list5$udata,"permat"= list5$permat, "persim" = list5$persim, "br" = br, "cl" = cl, "met" = met, "ep"= ep, "clep" = clep,"nn" = list5$nn, "sumper" = list5$sumper, "sumper2" = list5$sumper2, "ela" = list5$ela, "cel" = list5$cel,"endper" = list5$endper, "last" = list5$last)
  Matrices(result5,FALSE,let,sim,d,algo,algocol)
}


# A function that generates a geometric sequence
geomSeq <- function(start,ratio,begin,end){
  begin=begin-1
  end=end-1
  start*ratio**(begin:end)
}


Den <- function(lev,df,lastlist,flagtic){
  if(flagtic == TRUE) tic(sprintf("Den -- Level: %d ", lev))
  df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
  df_args <- c(df[str_which(names(df), "level.")], sep="/")
  if(lev == max(lastlist$clep)){
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
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  plot(xN)
}

# Create custom colour scheme
cs1 = make_col_scheme(chars=c("F","W","A","I","L","V","M","C","P","G","Y","T","S","H","K","R","E","D","Q","N"),
                      cols=c("#1E90FF", "#BA55D3", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#C6E2FF", "#C6E2FF", "#FFD700", "#00EE00", "#C1FFC1", "#54FF9F", "#54FF9F", "#FF0000", "#FF0000", "#FF0000", "#FFD700", "#FFD700", "#ED9121", "#ED9121"))

# A function to plot with logo level lev
LogoLev <- function(lev,lastlist,df,flagtic){
  if(flagtic == TRUE) tic(sprintf("LogoLev -- Level: %d ", lev))
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
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  ggseqlogo(listff, ncol=nc, method = "prob",col_scheme=cs1)
}

# Creating a plot for cluster 5 for example
LogoCl <- function(cl,lastlist,df,flagtic){
  if(flagtic == TRUE) tic(sprintf("LogoCl -- Cluster: %d ", cl))
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  ggseqlogo(na.omit(df[df[names(df) == sprintf('level.%d', lastlist$clep[cl])] == cl,]$AA.JUNCTION), method = "prob", col_scheme=cs1)
}

# A function to plot with barplot level lev
BarLev <- function(lev,lastlist,perlist,persimlist,Clus,let,sim,cho,flagtic){
  if(flagtic == TRUE) tic(sprintf("BarLev -- Level: %d ", lev))
  if(cho == "Identity"){
    t2 = which(lastlist$clep == lev)
    if(length(t2) %% 3 == 0){
      par(mfrow = c(length(t2)%/%3,3))
    }else{
      par(mfrow = c(length(t2)%/%3 + 1,3))
    }
    for(i in 1:length(t2)){
      ar = str_which(names(perlist),as.character(t2[i]))[1]
      par(xpd=TRUE)
      output <- matrix(unlist(perlist[ar]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
      barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d', t2[i])) 
      legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=let,cex = 0.6)
    }
  }else{
    t2 = which(lastlist$clep == lev)
    if(length(t2) %% 3 == 0){
      par(mfrow = c(length(t2)%/%3,3))
    }else{
      par(mfrow = c(length(t2)%/%3 + 1,3))
    }
    for(i in 1:length(t2)){
      ar = str_which(names(persimlist),as.character(t2[i]))[1]
      par(xpd=TRUE)
      output <- matrix(unlist(persimlist[ar]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
      barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d', t2[i])) 
      legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=names(sim),cex = 0.6)
    }
  }
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
}


# An alternative plot for cluster 3 (We must determine the cluster's permat)
BarCl <- function(cl,perlist,persimlist,Clus,let,sim,cho,lastlist,flagtic){
  if(flagtic == TRUE) tic(sprintf("BarCl -- Cluster :%d ", cl))
  if(cho == "Identity"){
    par(mar=c(3,3,4,4),xpd=TRUE)
    output <- matrix(unlist(perlist[cl]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
    barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d',cl)) 
    legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=let,cex = 0.6)
  }else{
    par(mar=c(3,3,4,4),xpd=TRUE)
    output <- matrix(unlist(persimlist[cl]), ncol = str_length(lastlist$udata$AA.JUNCTION[1]), byrow = FALSE)
    barplot(output[-nrow(output),], col=heat.colors(length(output[,1])-1), width=2, main = sprintf('Cluster.%d',cl)) 
    legend("topright",inset=c(-0.03,0), fill=heat.colors(length(output[,1])-1), legend=names(sim),cex = 0.6)
  }
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
}

# Sequences and Id's for specific level
AminoLev <- function(level,lastlist,df,Clus,flagtic){
  if(flagtic == TRUE) tic(sprintf("AminoLev -- Level: %d ", level))
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

# Sequences and Id's for specific cluster
AminoCl <- function(clust,lastlist,df,flagtic){
  #if(flagtic == TRUE) tic("AminoCl")
  if(clust == 0){
    x3 = data.frame("Sequence.ID" = df$Sequence.ID, "AA.JUNCTION" = df$AA.JUNCTION)
  }else{
    #lastlist$clep[clust]
    x3 = data.frame(na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$Sequence.ID), na.omit(df[df[which(names(df) == sprintf("level.%d", lastlist$clep[clust]))] == clust, ]$AA.JUNCTION))
    names(x3) = c("Sequence.ID","AA.JUNCTION")
    se = x3$Sequence.ID
    aa = as.character(x3$AA.JUNCTION)
    aa = as.list(aa)
    names(aa) = se 
  }
  #if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  x3 # Print in console
}


# A function which visualize the sequences of a data frame using common letters (i.e. "A _ _ _ _ K R _ _ _ Q _ Y Y Y _ _ _ T _")
Opt <- function(df,flagtic){
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
  paste(xar,collapse = ' ')
}

OptNew <- function(df){
  xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
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
    }else{
      xar[i] = names(which(f == max(f)))
    }
  }
  paste(xar,collapse = ' ')
}

# A function which visualize the sequences of a data frame using similarity groups (i.e. "A _ _ _ _ K Am _ Ba _ Ac _ Y Y Y _ _ _ T _")
Opt2 <- function(df,sim,alt,altsim,flagtic){
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

Opt2New <- function(df,sim){
  xar <- matrix(0,nrow=1, ncol=str_length(df$AA.JUNCTION[1]))
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
    }else{
      xar[i] = names(which(f == max(f)))
    }
  }
  paste(xar,collapse = ' ')
}

Id <- function(ff,cho,flagtic){
  if(flagtic == TRUE) tic(sprintf("Id -- Type: $s ", cho))
  if(cho == "Identity"){
    par(xpd=TRUE)
    pp <- data.frame(x = ff$level,y = ff$sumper)
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    ggplot(pp,aes(x = x,y = y)) + stat_sum()
  }else{
    par(xpd=TRUE)
    pp <- data.frame(x = ff$level,y = ff$sumper2)
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    ggplot(pp,aes(x = x,y = y)) + stat_sum()
  }
}


cc <- function(df,Clus,flagtic){
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
    fill = c("jj",ii),
    width = 1820,
    height = 775,
    collapsed = FALSE
  )
}

idenlev <- function(lev,Clus,cho,flagtic){
  if(flagtic == TRUE) tic(sprintf("idenlev -- Level: %d ", lev))
  if(cho == "Identity"){
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    sum(Clus[Clus$level == lev,]$Identity) /  length(Clus[Clus$level == lev,]$Identity)
    
  }else{
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    sum(Clus[Clus$level == lev,]$Similarity) /  length(Clus[Clus$level == lev,]$Similarity)
  }
}

idencl <- function(cl,Clus,cho,flagtic){
  if(flagtic == TRUE) tic(sprintf("idencl -- Cluster :%d ", cl))
  if(cho == "Identity"){
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    Clus[Clus$ClusterId == cl,]$Identity
  }else{
    if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
    Clus[Clus$ClusterId == cl,]$Similarity 
  }
}

EmPin <- function(lastlist,Clus,dfsd,flagtic){
  if(flagtic == TRUE) tic("EmPin")
  for (i in 0:(max(lastlist$clep))) {
    t1 = sum(Clus[Clus$level == i,]$Identity) /  length(Clus[Clus$level == i,]$Identity)
    t2 = sd(Clus[Clus$level == i,]$Identity)
    t3 = sum(Clus[Clus$level == i,]$Similarity) /  length(Clus[Clus$level == i,]$Similarity)
    t4 = sd(Clus[Clus$level == i,]$Similarity)
    rows = sprintf("level.%d", i)
    dfsd[i+1,] = c(rows,t1,t2,t3,t4)
  }
  colnames(dfsd)[1]="Level"
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  dfsd
}

Netw <- function(lev,thr,thrt,netyp,df,lastlist,Clus,sim,altsim,ts,shth,net_sil,flagtic){
  if(flagtic == TRUE) tic(sprintf("Network -- level: %d, Network type: %s, Threshold type: %s, Threshold value: %d, Using silhouette: %d ", lev,netyp,thrt,thr,net_sil))
  # We change the plotting parameters
  df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
  df_args <- c(df[str_which(names(df), "level.")], sep="/")
  if(lev == max(lastlist$clep)){
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
  ffg2sim = matrix(NA,nrow = max(bo)+1,ncol = (max(bo)+1))
  ffg3 = matrix(NA,nrow = max(bo)+1,ncol = (max(bo)+1)) # the final distance combining ffg and ffg2
  
  # anw trigwnikos
  for(i in 1:(max(bo)+1)){
    tempN1 = FindNode(xN,(sprintf("%d", i-1)))
    tempN1 = tempN1$path
    temp2N1 = Opt(AminoCl(i-1,lastlist,df,flagtic),flagtic)
    temp3N1 = Opt2(AminoCl(i-1,lastlist,df,flagtic),sim,TRUE,altsim,flagtic)
    for(j in i:(max(bo)+1)){          # max(df$clusters)+1
      tempN2 = FindNode(xN,(sprintf("%d", j-1)))
      tempN2 = tempN2$path
      tem = which(tempN1 == tempN2)
      tempN1[tem[length(tem)]]
      d1 = length(tempN1) - tem[length(tem)]
      d2 = length(tempN2) - tem[length(tem)]
      ffg[i,j]= d1 + d2
      temp2N2 = Opt(AminoCl(j-1,lastlist,df,flagtic),flagtic)
      temp2N2 = str_replace_all(temp2N2,"_"," ")
      temp3N2 = Opt2(AminoCl(j-1,lastlist,df,flagtic),sim,TRUE,altsim,flagtic)
      temp3N2 = str_replace_all(temp3N2,"_"," ")
      ffg2[i,j] = stringdist(temp2N1,temp2N2,method = "lv" )
      ffg2sim[i,j] = stringdist(temp3N1,temp3N2,method = "lv" )
      if( i == j ){
        ffg2[i,j] = str_length(lastlist$udata$AA.JUNCTION[1]) #gia na mhn fainontai velakia pisw
        ffg2sim[i,j] = str_length(lastlist$udata$AA.JUNCTION[1])
      }
    }
  }
  
  #find max
  ma =max(na.omit(as.vector(ffg)))
  
  ffg[which(ffg == 0)] = ma #gia na mhn fainontai velakia pisw
  ffg2 = str_length(lastlist$udata$AA.JUNCTION[1]) - ffg2
  ffg2sim = str_length(lastlist$udata$AA.JUNCTION[1]) - ffg2sim
  ffg3 = (2/str_length(lastlist$udata$AA.JUNCTION[1])) * ffg2 + (2/str_length(lastlist$udata$AA.JUNCTION[1])) * ffg2sim + (2 / ma) * (ma - ffg)
  
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
  tempor[which(thrtyp > thr)] = 0
  tempor2[which(thrtyp > thr)] = 0
  tempor2sim[which(thrtyp > thr)] = 0
  tempor3[which(thrtyp> thr)] = 0
  
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
    jhj = silhouette(Clus$level[1:(max(bo)+1)],t(tempor))
    matches <- regmatches(unique(x$pathString), gregexpr("[[:digit:]]+", unique(x$pathString)))
    for(i in 1:length(unique(x$pathString))){
      ttt =jhj[as.numeric(unlist(matches[i])),3]
      flag = FALSE
      j = length(ttt)
      mcl = length(ttt)
      while(flag == FALSE){
        if(ttt[j] >= ttt[j-1]){
          tempor3[j-1,j] = 0
          ffg3[j-1,j] = 0
          thrtyp[j-1,j] = 0
          thrtyp2[j-1,j] = 0
          j = j - 1
        }else{
          mcl = j
          flag = TRUE
        }
        
        if(j == 1){
          flag = TRUE
        }
      }
    }
  }
  
  jjj = matrix(1,nrow = max(bo)+1,ncol = (max(bo)+1))
  net0 <- graph_from_adjacency_matrix(jjj,mode = "upper")
  net1 <- graph_from_adjacency_matrix(jjj,mode = "upper")
  
  bb = vector(length = (max(bo)+1))
  bb[1]=0
  bb[2:length(bb)] = lastlist$clep[1:max(bo)]
  #XRWMA
  # Generate colors based on media type:
  colrs <- c("#1E90FF", "#BA55D3", "#0000FF", "#557fd2", "#54d17e", "#8aad62", "#C6E2FF", "#e5e234", "#FFD700", "#00EE00", "#C1FFC1", "#ea8509", "#54FF9F", "#FF0000", "#ed3b1c", "#ed1c7a", "#0c0c0c", "#b8d8af", "#ED9121","#45f713")
  colrs = colrs[1:(max(bb)+1)]
  V(net0)$color <- colrs[bb+1]
  V(net1)$color <- colrs[bb+1]
  
  # Compute node degrees (#links) and use that to set node size:
  deg <- (Clus$seqnum[1:length(bb)] / Clus$seqnum[1]) * 20
  V(net0)$size <- deg
  V(net1)$size <- deg
  
  # The labels are currently node IDs.
  # Setting them to NA will render no labels:
  V(net0)$label <- V(net0)-1
  V(net1)$label <- V(net1)-1
  
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
    
  net1.copy <- igraph::delete.edges(net1, which(E(net1)$width == 0))
    
  if(flagtic == TRUE) toc(log = TRUE,quiet = TRUE)
  listnet = list("net0.copy" = net0.copy, "net1.copy" = net1.copy, "bb" = bb, "colrs" = colrs, "ma" = ma, "ffg" = ffg, "ffg2" = ffg2, "ffg2sim" = ffg2sim, "ffg3" = ffg3 ,"tempor" = tempor, "tempor2" = tempor2, "tempor2sim" = tempor2sim, "tempor3" = tempor3,"thrtyp" = thrtyp, "thrtyp2" = thrtyp2)
  return(listnet)
}
# Insert Data
library("stringr")
library("dplyr")
# data <- read.csv(file.choose(), header = TRUE, sep = ";")
data <- read.csv("data/SampleData.csv", header = TRUE, sep = ";")


# Usefull Data
udata <- data[c(1,3)]
udata$AA.JUNCTION <- as.character(udata$AA.JUNCTION)

# Convert character to ascii
asc <- function(x) { strtoi(charToRaw(x),16L) }
asc("a")

# Convert ascii to character
chr <- function(n) { rawToChar(as.raw(n)) }
chr(asc("a"))

MyMatrix <- function(z,y){ # A function which computes the absolute and percentage matrix for a dataframe
  mymat <- matrix(0,nrow=20, ncol=26)
  permat <- matrix(0,nrow=20, ncol=26)
  
  gram <- 65          # letter "?" in ascii code
  for (i in 1:26){
    for (j in 1:20){
      if(i>1 && i<26 && rr != 0){
        permat[j,i-1] = (mymat[j,i-1] / rr) * 100 # percentege matrix
      }
      for (k in 1:length(z$Sequence.ID)){
        xar <- str_sub(z$AA.JUNCTION[k],j,j)
        if ( asc(xar) == gram){
          mymat[j,i] <- mymat[j,i] + 1 # absolute matrix
        }
      }
    }
    rr = sum(mymat[,i]) # sum for every row
    gram <- gram + 1
  }
  # completing the percentage matrix with the last letter 
  if(rr != 0){
    for(i in 1:20){
      permat[i,26] = (mymat[i,26] / rr) * 100
    }
  }
  # Deciding which matrix we want to return
  if (y == 0){
  return (mymat)
  }
  return (permat)
}

# Create the absolute matrix for the udata
umat = MyMatrix(udata,0)

# Create the percentege matrix for the udata
uper = MyMatrix(udata,1)

# finding the 2 sub-data frames (example for cell(11,18))
l = chr(64+18) # letter
x = filter(udata,str_detect(str_sub(udata$AA.JUNCTION,11,11), l))
y = filter(udata,!str_detect(str_sub(udata$AA.JUNCTION,11,11), l))

# finding the percentage matrix for these data frames 
per1 = MyMatrix(x,1)
per2 = MyMatrix(y,1)

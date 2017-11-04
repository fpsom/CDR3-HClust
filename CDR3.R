# Insert Data
data <- read.csv(file.choose(), header = TRUE, sep = ";")

# Usefull Data
udata <- data[c(1,3)]
udata$AA.JUNCTION <- as.character(udata$AA.JUNCTION)

# Convert character to ascii
asc <- function(x) { strtoi(charToRaw(x),16L) }
asc("a")

# Create a null matrix(20,26)
mymat <- matrix(0,nrow=20, ncol=26)

# Create the absolute matrix
for (i in 1:20){
  gram <- 65          # letter "Á" in ascii code
  for (j in 1:26){
    for (k in 1:length(udata$Sequence.ID)){
      xar <- strsplit(udata$AA.JUNCTION[k], '')[[1]][i]
      if ( asc(xar) == gram){
        mymat[i,j] <- mymat[i,j] + 1
      }
    }
    gram <- gram + 1
  }
}
mymat

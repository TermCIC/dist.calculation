# Calculating dna p-dsitances and tree distances
# Progammer: Chun-I Chiu
install.packages(reshape2)
install.packages(ape)
library(ape)
library(reshape2)

dist.calculation <-
  function(treeFileName,
           seqFileName,
           savingFileName) {
    # Reading the dna sequence & calculating p-distance
    seqs <- read.nexus.data(seqFileName) # Your input
    write.dna(seqs,
              file = "temp.fas",
              format = "fasta",
              append = FALSE)
    seqs <- read.dna("temp.fas", format = "fasta")
    pDist <-
      dist.dna(
        seqs,
        model = "raw",
        as.matrix = T,
        pairwise.deletion = T
      )
    cNames <- colnames(pDist)
    rNames <- rownames(pDist)
    pDist <- cbind(cNames, pDist)
    pDist <- pDist[order(pDist[, 1]),]
    pDist <- pDist[,-1]
    pDist <- rbind(rNames, pDist)
    pDist <- pDist[,order(pDist[1, ])]
    pDist <- pDist[-1,]
    taxonLength <- nrow(pDist)
    counter1 <- 1
    counter2 <- 1
    for (j in 1:taxonLength) {
      for (k in 1:taxonLength) {
        if (k <= counter2) {
          pDist[j, k] <- NA
        } else {
          pDist[j, k] <- pDist[j, k]
          counter1 <- counter1 + 1
        }
      }
      counter2 <- counter2 + 1
    }
    pDist.col <- melt(pDist)
    pDist.col <- pDist.col[complete.cases(pDist.col), ]
    pDist.col <- pDist.col[order(pDist.col[, 1], pDist.col[, 2]), ]
    # Reading the tree & calculating distances
    tree <- read.nexus(treeFileName)
    MBDist <- cophenetic(tree)
    cNames <- colnames(MBDist)
    rNames <- rownames(MBDist)
    MBDist <- cbind(cNames, MBDist)
    MBDist <- MBDist[order(MBDist[, 1]),]
    MBDist <- MBDist[,-1]
    MBDist <- rbind(rNames, MBDist)
    MBDist <- MBDist[,order(MBDist[1, ])]
    MBDist <- MBDist[-1,]
    counter1 <- 1
    counter2 <- 1
    for (j in 1:taxonLength) {
      for (k in 1:taxonLength) {
        if (k <= counter2) {
          MBDist[j, k] <- NA
        } else {
          MBDist[j, k] <- MBDist[j, k]
          counter1 <- counter1 + 1
        }
      }
      counter2 <- counter2 + 1
    }
    MBDist.col <- melt(MBDist)
    MBDist.col <- MBDist.col[complete.cases(MBDist.col), ]
    MBDist.col <- MBDist.col[order(MBDist.col[, 1], MBDist.col[, 2]), ]
    Table.length <- nrow(MBDist.col)
    Output.table <- cbind(MBDist.col, pDist.col)
    Output.table <- Output.table[,-(4:5)]
    colnames(Output.table) <-
      c("Species 1", "Species 2", "Tree distance", "DNA p-distance")
    write.csv(Output.table, savingFileName)
  }
#------------------------------------------------------------------------------------------
dist.calculation("Documents/Global_Dryocosmus/2020/saturation_plotting/concatenated.nex", # Input tree
                 "Documents/Global_Dryocosmus/2020/saturation_plotting/pos3.nexus", # Input sequence
                 "Documents/Global_Dryocosmus/2020/saturation_plotting/distTable.csv") # Output the csv file


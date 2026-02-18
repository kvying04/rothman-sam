library(ggplot2)

#SET PATH

# Set your working directory to the Desktop
setwd("~/Desktop/slidingwindow")
getwd() 



#path <- file.path("/Users/samanthafiallo/Downloads/slidingwindow")
#setwd(path)


#SET TERMS
chromosomes <- c("I", "II", "III", "IV", "V", "X", "mtDNA") #setting constants
cLengths <- c(16000000, 16000000, 14000000, 18000000, 21000000, 18000000, 14000)
cxGroups <- c("JR4295", "JR4296", "JR4297") #MANIPULATE INDV INCLUDED HERE

#CREATE DATA FRAMES
#chromosomes_filenames <- c("chrI.csv", "chrII.csv", "chrIII.csv", "chrIV.csv", "chrV.csv", "chrX.csv") DON'T THINK NEED FOR 
#chromosomes_filenames <- c("chrI.csv") #TO TRY ONE AT A TIME
chrAll <- read.csv("windowGeneration/withreference_4_reformatted_forslidingWindow.csv", sep = ',' )    
for(c in 1:7){
  chrName <- chromosomes[c]
  print(chrName)
  chrLen <- cLengths[c]
  #chrAll <- read.csv("windowGeneration/withreference_4_reformatted_forslidingWindow.csv", sep = ',' ) #CHANGE CSV HERE
  chr<-subset(chrAll, CHROM==chrName) 
  print("howdy")
  
  #this sets the window 
  window = as.integer(1)
  
  position=c()
  snpCountcx=c()
  cxCount <- 0
  #for(p in seq(5820000, 5850000, by=1)) {  # adjustable 
  for(p in seq(window, chrLen, by=100000)) {  # Endpoint should be adjusted
    #for(p in seq(window, 5000, by=1000)) {  # test case
    lbound <- p - window
    hbound <- p + window
    cxCount <- 0
    
    # Filter data within bounds
    inWindow <- subset(chr, CHROM >= lbound & CHROM <= hbound)
    
    # Loop through SNPs in the window
    for(snp in rownames(inWindow)) {
      
      # Check for low group SNPs
      # Check for high group SNPs
      for(highg in highGroups) {
        if(substring(inWindow[snp, highg], 1, 3) == "1/1") {
          hCount <- hCount + 1
        }
      }
    }
    
    # Update the vectors with the counts for this window
    # snpCountL <- append(snpCountL, lCount)
    # snpCountH <- append(snpCountH, hCount)
    # snpCountD <- append(snpCountD, lCount-hCount)
    # position <- append(position, p)
    
    # Update the vectors with the counts for this window
    # snpCountL <- append(snpCountL, lCount/length(lowGroups))
    # snpCountH <- append(snpCountH, hCount/length(highGroups))
    # snpCountD <- append(snpCountD, lCount-hCount)
    # position <- append(position, p)
    
    snpCountL <- c(snpCountL, lCount/length(lowGroups))
    snpCountH <- c(snpCountH, hCount/length(highGroups))
    snpCountD <- c(snpCountD, lCount - hCount)
    position <- c(position, p)
  }
  
  hist_data <- data.frame(position, snpCountL, snpCountH, snpCountD)
 #SUBTRACTION GRAPH
   ggplot(data=hist_data) +
    geom_line(aes(x=position, y=snpCountD), color="black") +
    
    theme_classic() +
    ggtitle(chrName)+
    ylab("SNP Count") +
    xlab("Position") +
    
    theme(
      axis.title.x = element_text(size=20),
      axis.title.y = element_text(size=20),
      axis.text.x = element_text(size=20),
      axis.text.y = element_text(size=20),
      legend.text = element_text(size=20),
      legend.title = element_text(size=20)
    )
  print(last_plot())
}

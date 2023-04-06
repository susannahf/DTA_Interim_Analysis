# create test data for use in analysis
# S Fleming 2022
source("generateDTAdata.R")


# create a 2x2 table with Se=65%, Sp=85%, prevalence=35%, and 1000 samples
tbt <- create2x2(n=1000, p=0.35, se=0.65, sp=0.85, digits=4, verbose=TRUE)


# the following code creates two example datasets from the 2x2 table
# the data points are randomised to provide different example "studies"
# setting the random seed ensures that the same dataset is created each time
# the code is run

set.seed(2022) # Se and Sp both start extreme and settle normally
testData1r2022 <- DTAdataFrom2x2(tbt, randomise = T)

# add continuous sensitivity and specificity to data
contsespr <- continuousSeSp(testData1r2022)

#plot
x=1:1000
plot(x,contsespr$Se,'l',col="red",ylim=c(0,1), main="Seed of 2022, testData1")
lines(x,contsespr$Sp,col="blue")
  

# second dataset
set.seed(9999) # Se and Sp both start very high, drop, peak, then settle
testData1r9999 <- DTAdataFrom2x2(tbt, randomise = T)
contsespr <- continuousSeSp(testData1r9999)

contsespr

#plot
x=1:1000
plot(x,contsespr$Se,'l',col="red",ylim=c(0,1), main="Seed of 9999, testData2")
lines(x,contsespr$Sp,col="blue")

# save objects to files for later use
saveRDS(testData1r2022, "testData1.rds")
saveRDS(testData1r9999, "testData2.rds")
# example code for how to run DTA interim analysis

source("generateDTAdata.R") # for continuousSeSp()
source("DTAinterimAnalysis.R") # for DTAdiscreteInterimAnalysis


# create a (not completely) random data set to practice on
testdata <- data.frame( 
  reference = as.logical(rbinom(100,1,0.2)),
  index = as.logical(rbinom(100,1,0.4))
)

# continuousSeSp adds TP and TN columns, which are needed for interim analysis code
testwithTPTN <- continuousSeSp(testdata)

# run interim analyses at n=20, 50, and 100
# assuming that the planned sample size is 150
# and cutoffs of 0.6 for both Se and Sp
# the minimal set of columns for this code are "reference", "TP", and "TN"
print("Interim analysis at 20, 50, and 100 points")
test1 <- DTAdiscreteInterimAnalysis(testwithTPTN,c(20,50,100), pSe=0.6, pSp=0.6, N=150, prevalence = 0.2)
print(test1)

print("Continuous interim analysis up to 100 points")
test2 <- DTAdiscreteInterimAnalysis(testwithTPTN,1:100, pSe=0.6, pSp=0.6, N=150, prevalence = 0.2)
print(test2)

# you can also run interim analyses using just the cumulative counts
# this is the data from dataset2
print("Interim analysis using cumulative counts")
cumulativedata <- data.frame(N=c(20,50,100),
                             RefT=c(13,29,51),
                             TP=c(9, 22, 35),
                             TN=c(6, 19, 51))
test3 <- DTAcumulativeInterimAnalysis(cumulativedata, pSe=0.8, pSp=0.9, N=200, prevalence = 0.4)
print(test3)







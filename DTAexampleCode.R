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
# and cutoffs of 0.6 for both Se and Sp
# the minimal set of columns for this code are "reference", "TP", and "TN"
test1 <- DTAdiscreteInterimAnalysis(testwithTPTN,c(20,50,100), pSe=0.6, pSp=0.6, N=100, prevalence = 0.2, simpleOutput = FALSE)
print(test1)

test2 <- DTAdiscreteInterimAnalysis(testwithTPTN,1:100, pSe=0.6, pSp=0.6, N=100, prevalence = 0.2, simpleOutput = FALSE)
print(test2)



# Trying to debug DTAInterimAnalysis
source("generateDTAdata.R") # for continuousSeSp()
source("DTAinterimAnalysis.R") # for DTAdiscreteInterimAnalysis

# let's try one of the sample data sets
testdata1 <- readRDS("testData1.rds")

# continuousSeSp adds TP and TN columns, which are needed for interim analysis code and shortens to 100pts
testdata1 <- continuousSeSp(testdata[1:100,])

# run interim analyses at n=20, 50, and 100
# and cutoffs of 0.6 for both Se and Sp
# the minimal set of columns for this code are "reference", "TP", and "TN"
test1 <- DTAdiscreteInterimAnalysis(testwithTPTN,c(20,50,100), pSe=0.6, pSp=0.6, simpleOutput = TRUE)
print(test1)

test2 <- DTAdiscreteInterimAnalysis(testwithTPTN,1:100, pSe=0.6, pSp=0.6, simpleOutput = TRUE)
print(test2)


# debugging
termthresh <- modelFlemingTerminationThresholds(0.4, n=100)

#Se
plot(x,testdata1$Se,'l',col="red",ylim=c(0,1), main="Continuous Se", ylab="Sensitivity")
lines(x, termthresh$invH0limit, col="red", 'l', lty="dashed")
lines(x, termthresh$invH1limit, col="red", 'l', lty="dotdash")
lines(x, rep(0.6,100), col="red")
grid()
#legend("bottomright", c("Sensitivity", "H0 limit", "H1limit", "p0"), col=c("red"),lty=c(1,2,4,1))
#Sp
plot(x,testdata1$Sp,'l',col="blue",ylim=c(0,1), main="Continuous Sp", ylab="Specificity")
lines(x, termthresh$invH0limit, col="blue", 'l', lty="dashed")
lines(x, termthresh$invH1limit, col="blue", 'l', lty="dotdash")
lines(x, rep(0.6,100), col="blue")
grid()
#legend("bottomright", c("Specificity", "H0 limit", "H1limit", "p0"), col=c("blue"),lty=c(1,2,4,1))




















# # create a (not completely) random data set to practice on
# set.seed(7)
# testdata <- data.frame( 
#   reference = as.logical(rbinom(100,1,0.2)),
#   index = as.logical(rbinom(100,1,0.4))
# )
# 
# # continuousSeSp adds TP and TN columns, which are needed for interim analysis code
# testwithTPTN <- continuousSeSp(testdata)
# 
# # run interim analyses at n=20, 50, and 100
# # and cutoffs of 0.6 for both Se and Sp
# # the minimal set of columns for this code are "reference", "TP", and "TN"
# test1 <- DTAdiscreteInterimAnalysis(testwithTPTN,c(20,50,100), pSe=0.6, pSp=0.6, simpleOutput = TRUE)
# print(test1)
# 
# test2 <- DTAdiscreteInterimAnalysis(testwithTPTN,1:100, pSe=0.6, pSp=0.6, simpleOutput = TRUE)
# print(test2)
# 
# 
# # debugging
# termthresh <- modelFlemingTerminationThresholds(0.4, n=100)
# 
# #Se
# plot(x,testwithTPTN$Se,'l',col="red",ylim=c(0,1), main="Continuous Se", ylab="Sensitivity")
# lines(x, termthresh$invH0limit, col="red", 'l', lty="dashed")
# lines(x, termthresh$invH1limit, col="red", 'l', lty="dotdash")
# lines(x, rep(0.6,100), col="red")
# grid()
# #legend("bottomright", c("Sensitivity", "H0 limit", "H1limit", "p0"), col=c("red"),lty=c(1,2,4,1))
# #Sp
# plot(x,testwithTPTN$Sp,'l',col="blue",ylim=c(0,1), main="Continuous Sp", ylab="Specificity")
# lines(x, termthresh$invH0limit, col="blue", 'l', lty="dashed")
# lines(x, termthresh$invH1limit, col="blue", 'l', lty="dotdash")
# lines(x, rep(0.6,100), col="blue")
# grid()
# #legend("bottomright", c("Specificity", "H0 limit", "H1limit", "p0"), col=c("blue"),lty=c(1,2,4,1))

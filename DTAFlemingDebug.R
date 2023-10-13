# something strange is happening, where sensitivity in particular is terminating with the "wrong" conclusion.

# but it seems to work on the test data in teh RMD file

# plan: 
# replicate the test data @ 200 points results which seem to be correct
# try running the test data @ 100 points

source("DTAinterimAnalysis.R")
source("generateDTAdata.R")

# load test data
testdata1 <- readRDS("testData1.rds")
testdata2 <- readRDS("testData2.rds")

# set up short test data
testshort1 <- testdata1[1:200,]
testshort2 <- testdata2[1:200,]
# these should nominally have Se= 65 and Sp = 85 but won't 

# calculate continuous Se and Sp
contshort1 <- continuousSeSp(testshort1)
contshort2 <- continuousSeSp(testshort2)
contshorter1 <- contshort1[1:100, ]
contshorter2 <- contshort2[1:100, ]

# find the theoretical Fleming thresholds for these
sep0_60 <- modelFlemingTerminationThresholds(0.4, n=200)
spp0_90 <- modelFlemingTerminationThresholds(0.1, n=200)
sep0_60short <- modelFlemingTerminationThresholds(0.4, n=100)
spp0_90short <- modelFlemingTerminationThresholds(0.1, n=100)

# plot data against thresholds
x <- 1:nrow(contshort1)
xs <- 1:nrow(contshorter1)

# plot continuous Se and Sp for contshort1
#Se 200
plot(x,contshort1$Se,'l',col="red",ylim=c(0,1), main="Continuous Se, data set 1, n=200, p0=0.6", ylab="Sensitivity")
lines(x, sep0_60$invH0limit, col="red", 'l', lty="dashed")
lines(x, sep0_60$invH1limit, col="red", 'l', lty="dotdash")
lines(x, rep(0.6,200), col="red")
grid()
#Se 100
plot(xs,contshorter1$Se,'l',col="red",ylim=c(0,1), main="Continuous Se, data set 1, n=100, p0=0.6", ylab="Sensitivity")
lines(xs, sep0_60short$invH0limit, col="red", 'l', lty="dashed")
lines(xs, sep0_60short$invH1limit, col="red", 'l', lty="dotdash")
lines(xs, rep(0.6,100), col="red")
grid()
#Sp 200
plot(x,contshort1$Sp,'l',col="blue",ylim=c(0,1), main="Continuous Sp, data set 1, n=200, p0=0.9", ylab="Specificity")
lines(x, spp0_90$invH0limit, col="blue", 'l', lty="dashed")
lines(x, spp0_90$invH1limit, col="blue", 'l', lty="dotdash")
lines(x, rep(0.9,200), col="blue")
grid()
#Sp 100
plot(xs,contshorter1$Sp,'l',col="blue",ylim=c(0,1), main="Continuous Sp, data set 1, n=100, p0=0.9", ylab="Specificity")
lines(xs, spp0_90short$invH0limit, col="blue", 'l', lty="dashed")
lines(xs, spp0_90short$invH1limit, col="blue", 'l', lty="dotdash")
lines(xs, rep(0.9,100), col="blue")
grid()
# plot continuous Se and Sp for contshort2
#Se 200
plot(x,contshort2$Se,'l',col="red",ylim=c(0,1), main="Continuous Se, data set 2, n=200, p0=0.6", ylab="Sensitivity")
lines(x, sep0_60$invH0limit, col="red", 'l', lty="dashed")
lines(x, sep0_60$invH1limit, col="red", 'l', lty="dotdash")
lines(x, rep(0.6,200), col="red")
grid()
#Se 100
plot(xs,contshorter2$Se,'l',col="red",ylim=c(0,1), main="Continuous Se, data set 2, n=100, p0=0.6", ylab="Sensitivity")
lines(xs, sep0_60short$invH0limit, col="red", 'l', lty="dashed")
lines(xs, sep0_60short$invH1limit, col="red", 'l', lty="dotdash")
lines(xs, rep(0.6,100), col="red")
grid()
#Sp 200
plot(x,contshort2$Sp,'l',col="blue",ylim=c(0,1), main="Continuous Sp, data set 2, n=200, p0=0.9", ylab="Specificity")
lines(x, spp0_90$invH0limit, col="blue", 'l', lty="dashed")
lines(x, spp0_90$invH1limit, col="blue", 'l', lty="dotdash")
lines(x, rep(0.9,200), col="blue")
grid()
#Sp 100
plot(xs,contshorter2$Sp,'l',col="blue",ylim=c(0,1), main="Continuous Sp, data set 2, n=100, p0=0.9", ylab="Specificity")
lines(xs, spp0_90short$invH0limit, col="blue", 'l', lty="dashed")
lines(xs, spp0_90short$invH1limit, col="blue", 'l', lty="dotdash")
lines(xs, rep(0.9,100), col="blue")
grid()


# now lets run continuous interim analysis to see if it terminates where and how we expect.
pSe=0.6
pSp=0.9
testpoints <- 1:200
test1 <- DTAdiscreteInterimAnalysis(contshort1,testpoints, pSe=pSe, pSp=pSp)
test2 <- DTAdiscreteInterimAnalysis(contshort2,testpoints, pSe=pSe, pSp=pSp)
testpoints <- 1:100
test1short <- DTAdiscreteInterimAnalysis(contshorter1,testpoints, pSe=pSe, pSp=pSp)
test2short <- DTAdiscreteInterimAnalysis(contshorter2,testpoints, pSe=pSe, pSp=pSp)

# report results
print("Test data 1, n=200:")
print(test1)
print("Test data 1, n=100:")
print(test1short)
print("Test data 2, n=200:")
print(test2)
print("Test data 2, n=100:")
print(test2short)

# this gives a few situations where it would appear actual stopping and predicted stopping do not match

#next plan
# 2) check the stopping point code
# 3) replicate the stopping using only the original Fleming code

# 1) visualise the plots around the stopping points
#Se 200
plot(x,contshort1$Se,'l',col="red",ylim=c(0,1), main="Continuous Se, data set 1, n=200, p0=0.6", ylab="Sensitivity",
     xlim=c(test1$Sensitivity$terminateAt-10, test1$Sensitivity$terminateAt+10))
lines(x, sep0_60$invH0limit, col="red", 'l', lty="dashed")
lines(x, sep0_60$invH1limit, col="red", 'l', lty="dotdash")
lines(x, rep(0.6,200), col="red")
grid(nx=2)
#Se 100
plot(xs,contshorter1$Se,'l',col="red",ylim=c(0,1), main="Continuous Se, data set 1, n=100, p0=0.6", ylab="Sensitivity",
     xlim=c(test1short$Sensitivity$terminateAt-10, test1short$Sensitivity$terminateAt+10))
lines(xs, sep0_60short$invH0limit, col="red", 'l', lty="dashed")
lines(xs, sep0_60short$invH1limit, col="red", 'l', lty="dotdash")
lines(xs, rep(0.6,100), col="red")
grid(nx=2)
#Sp 200
plot(x,contshort1$Sp,'l',col="blue",ylim=c(0,1), main="Continuous Sp, data set 1, n=200, p0=0.9", ylab="Specificity",
     xlim=c(test1$Specificity$terminateAt-10, test1$Specificity$terminateAt+10))
lines(x, spp0_90$invH0limit, col="blue", 'l', lty="dashed")
lines(x, spp0_90$invH1limit, col="blue", 'l', lty="dotdash")
lines(x, rep(0.9,200), col="blue")
grid(nx=2)
#Sp 100
plot(xs,contshorter1$Sp,'l',col="blue",ylim=c(0,1), main="Continuous Sp, data set 1, n=100, p0=0.9", ylab="Specificity",
     xlim=c(test1short$Specificity$terminateAt-10, test1short$Specificity$terminateAt+10))
lines(xs, spp0_90short$invH0limit, col="blue", 'l', lty="dashed")
lines(xs, spp0_90short$invH1limit, col="blue", 'l', lty="dotdash")
lines(xs, rep(0.9,100), col="blue")
grid(nx=2)
# plot continuous Se and Sp for contshort2
#Se 200
plot(x,contshort2$Se,'l',col="red",ylim=c(0,1), main="Continuous Se, data set 2, n=200, p0=0.6", ylab="Sensitivity",
     xlim=c(test2$Sensitivity$terminateAt-10, test2$Sensitivity$terminateAt+10))
lines(x, sep0_60$invH0limit, col="red", 'l', lty="dashed")
lines(x, sep0_60$invH1limit, col="red", 'l', lty="dotdash")
lines(x, rep(0.6,200), col="red")
grid(nx=2)
#Se 100
plot(xs,contshorter2$Se,'l',col="red",ylim=c(0,1), main="Continuous Se, data set 2, n=100, p0=0.6", ylab="Sensitivity",
     xlim=c(test2short$Sensitivity$terminateAt-10, test2short$Sensitivity$terminateAt+10))
lines(xs, sep0_60short$invH0limit, col="red", 'l', lty="dashed")
lines(xs, sep0_60short$invH1limit, col="red", 'l', lty="dotdash")
lines(xs, rep(0.6,100), col="red")
grid(nx=2)
#Sp 200
plot(x,contshort2$Sp,'l',col="blue",ylim=c(0,1), main="Continuous Sp, data set 2, n=200, p0=0.9", ylab="Specificity",
     xlim=c(test2$Specificity$terminateAt-10, test2$Specificity$terminateAt+10))
lines(x, spp0_90$invH0limit, col="blue", 'l', lty="dashed")
lines(x, spp0_90$invH1limit, col="blue", 'l', lty="dotdash")
lines(x, rep(0.9,200), col="blue")
grid(nx=2)
#Sp 100
plot(xs,contshorter2$Sp,'l',col="blue",ylim=c(0,1), main="Continuous Sp, data set 2, n=100, p0=0.9", ylab="Specificity",
     xlim=c(test2short$Specificity$terminateAt-10, test2short$Specificity$terminateAt+10))
lines(xs, spp0_90short$invH0limit, col="blue", 'l', lty="dashed")
lines(xs, spp0_90short$invH1limit, col="blue", 'l', lty="dotdash")
lines(xs, rep(0.9,100), col="blue")
grid(nx=2)

# stopping points do not correspond with the theory. They are close for the n=200, but not **identical**



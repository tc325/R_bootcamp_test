RandProjBinary1 <- function(p=1000,d=10)
    {
        R <- matrix(1/sqrt(p)*2*(rbinom(d*p,1,0.5)-0.5),p,d)
        R <- qr.Q(qr(R))[,1:d]
        return(R)
    }

RandProjBinary2 <- function(p=1000,d=10)
    {
        R <- matrix(sqrt(3)/sqrt(p)*sample(c(-1,0,1),p*d,replace = T,c(1/6,2/3,1/6)),p,d)
        R <- qr.Q(qr(R))[,1:d]
        return(R)
    }


RandProjGauss <- function(p=1000,d=10)
    {
        R <- matrix(1/sqrt(p)*rnorm(p*d,0,1),p,d)
        R <- qr.Q(qr(R))[,1:d]
        return(R)
    }



library(ascrda)
library(class)
library(mvtnorm)
library(MASS)
library(Matrix)
library(datamicroarray)
library(pamr)
library(penalizedLDA)


RPClassifierAlt <- function(XTrain, YTrain, XTest, YTest, d, J = 1000, Projection = "Gauss",M=5) 
{
    n = length(YTrain)
    n.test <- length(YTest)
    p = ncol(XTrain)
    n1 <- floor(n/5)
    Cutoff.valt1 <- matrix(0,M,n1)
    Cutoff.valt2 <- matrix(0,M,n1)
    Cutoff.val1 <- matrix(0,J,n1)
    Cutoff.val2 <- matrix(0,J,n1)
    Class.vote1 <- matrix(0,J,n.test)
    Class.vote2 <- matrix(0,J,n.test)
    S.n <- rep(0,J)
    kcv.vote <- rep(0,5)
    kcv.voteRP <- rep(0,5)
    k <- c(3,5,10,15,25)
    for (j in 1:J)
        {
            if (Projection=="Gauss") RP <- RandProjGauss(p,d)
            for (l in 1:5)
                {
                    kcv.voteRP[l] <= sum(knn.cv(XTrain%*%RP, YTrain, k[l]) != YTrain)
                }
            for (m in 1:M)
                {
                    s <- sample(1:n,n1)
                    knn.out <- knn(XTrain[-s,]%*%RP, XTrain[s,]%*%RP, YTrain[-s], k[order(kcv.voteRP)[1]],prob = TRUE)
                    Cutoff.valt2[m,] <- as.numeric(knn.out)
                }
            Cutoff.val2[j,] <- colMeans(Cutoff.valt2)
            knn.out2 <- knn(XTrain%*%RP, XTest%*%RP, YTrain,  k[order(kcv.voteRP)[1]], prob = TRUE)
            Class.vote2[j,] <- as.numeric(knn.out2)
        }
    Class2 <- 1 + as.numeric(colMeans(Class.vote2) > mean(Cutoff.val2))
    errRP2 <- 100*sum(Class2 != YTest)/n.test
    for (l in 1:5)
        {
            kcv.vote[l] <= sum(knn.cv(XTrain, YTrain, k[l]) != YTrain)
        }
    err <-  100*sum(as.numeric(knn(XTrain, XTest, YTrain, k[order(kcv.vote)[1]]) != YTest))/n.test
    return(list(Class2 = Class2, error2 = errRP2, error3 = err))
}






#data('golub', package = 'datamicroarray')


#data('sun', package = 'datamicroarray')

#Spam <- read.table("spambase.data",sep = ",")
#Ozone <- read.table("ozone.data",sep = ",")
#MiniBoo <- read.table("MiniBooNE_PID.txt")
#MiniBoo.Y <-c(rep(1, 36499),rep(2, 93565))

#bio.data <- read.table("biodeg.csv",sep = ";")

Model1 <- function(n,p,s0){
    Y1 <- rmultinom(1,n,c(1/2,1/2))
    Y <- c(rep(1,Y1[1,1]),rep(2,Y1[2,1]))
    mu1 <- c(rep(-1.56/sqrt(s0),s0),rep(0,p-s0))
    mu2 <- c(rep(1.56/sqrt(s0),s0),rep(0,p-s0))
    X1 <- mvrnorm(Y1[1,1],mu1,diag(p))
    X2 <- mvrnorm(Y1[2,1],mu2,diag(p))
    X  <- rbind(X1,X2)
    return(list(x=X,y=Y))
}

Model2 <- function(n,p,s0){
    Y1 <- rmultinom(1,n,rep(1/2,2))
    Y <- c(rep(1,Y1[1,1]),rep(2,Y1[2,1]))
    mu1 <- c(rep(1,s0),rep(0,p-s0))
    mu2 <- c(rep(1,s0),rep(0,p-s0))
    Sigma <- matrix(0.5,p,p) + 0.5*diag(1,p,p)
    mu1 <- -1.56*mu1/sqrt(t(mu1)%*%solve(Sigma)%*%mu1)
    mu2 <- 1.56*mu2/sqrt(t(mu2)%*%solve(Sigma)%*%mu2)
    X1 <- mvrnorm(Y1[1,1],mu1,Sigma)
    X2 <- mvrnorm(Y1[2,1],mu2,Sigma)
    X  <- rbind(X1,X2)
   return(list(x=X,y=Y))
}

L <- 100

errRP1 <- rep(0,100)
errRP2 <- rep(0,100)
errRP3 <- rep(0,100)
errknn <- rep(0,100)
errpenLDA <- rep(0,100)
errNSC <- rep(0,100)
errSCRDA <- rep(0,100)
errIR <- rep(0,100)

Model <- Model1(10000,50,1)

data.X <-  Model$x #dataset covariates 
data.Y <-  Model$y  #class labels
for (l in 1:L)
    {
        s <- sample(1:length(data.Y),1000,replace  = FALSE) #training set
        s1 <- sample((1:length(data.Y))[-s],1000,replace = FALSE) #test set
        out.RP5 <- RPClassifierAlt(data.X[s,],data.Y[s], data.X[s1,], data.Y[s1], d = 5, J = 1000, Projection = "Gauss", M = 5)
        out.RP10 <- RPClassifierAlt(data.X[s,],data.Y[s], data.X[s1,], data.Y[s1], d = 10, J = 1000, Projection = "Gauss", M = 5)
        out.RP20 <- RPClassifierAlt(data.X[s,],data.Y[s], data.X[s1,], data.Y[s1], d = 20, J = 1000, Projection = "Gauss", M = 5)

        errRP1[l] <- out.RP5$error2
        errRP2[l] <- out.RP10$error2
        errRP3[l] <- out.RP20$error2
      
        errknn[l] <- out.RP5$error3

        cv.penLDA <- PenalizedLDA.cv(data.X[s,],data.Y[s],lambdas=c(0.01,0.05,0.1,0.5,1))
        out.penLDA  <- PenalizedLDA(data.X[s,],data.Y[s],data.X[s1,],lambda = cv.penLDA$bestlambda , K = cv.penLDA$bestK)
        errpenLDA[l] <- 100*sum(out.penLDA$ypred[,1]-data.Y[s1] != 0)/length(data.Y[s1])

        train.NSC <- pamr.train(list(x = t(data.X[s,]), y = data.Y[s]))
        out.NSC <- pamr.predict(train.NSC, t(data.X[s1,]), threshold = 1, type = c("class"))
        errNSC[l] <- 100*sum(out.NSC != data.Y[s1])/length(data.Y[s1])

        out.SCRDA <- ascrda(data.X[s,], data.Y[s], data.X[s1,], data.Y[s1], SCRDAmethod = "SCRDA")
        errSCRDA[l] <- 100*as.numeric(out.SCRDA$SCRDA)

        errIR[l] <- FitDLDA(data.X[s,], data.Y[s], data.X[s1,], data.Y[s1])$Err*100
        
    }

mean(errRP1)
sqrt(var(errRP1))/sqrt(L)

mean(errRP2)
sqrt(var(errRP2))/sqrt(L)

mean(errRP3)
sqrt(var(errRP3))/sqrt(L)

mean(errknn)
sqrt(var(errknn))/sqrt(L)

mean(errpenLDA)
sqrt(var(errpenLDA))/sqrt(L)

mean(errNSC)
sqrt(var(errNSC))/sqrt(L)

mean(errSCRDA)
sqrt(var(errSCRDA))/sqrt(L)

mean(errIR)
sqrt(var(errIR))/sqrt(L)










#TestX <- read.table("madelon_test.data")
#TrainX <- read.table("madelon_train.data")
#ValX <- read.table("madelon_valid.data")
#Mad.X <- rbind(as.matrix(TrainX), as.matrix(ValX))
#TrainY <- read.table("madelon_train.labels")
#ValY <- read.table("madelon_valid.labels")
#Mad.Y <- c(1/2*(3-as.numeric(TrainY[,1])),1/2*(3-as.numeric(ValY[,1])))


#Nomao <- read.csv("Nomao.data",sep = ",")
#Nomao.X <- data.matrix(Nomao[,c(2:7,10:15,18:23,26:31,34:39,42:47,50:55,58:63,66:71,74:79,82:87,90:92,94:96,98:100,102:104,106:108,110:112,114:116,118:119)])
#Nomao.Y <- 1/2*(3-as.numeric(Nomao[,120]))
#Nomao.X <- Nomao.X[order(Nomao.Y),]
#Nomao.Y <- Nomao.Y[order(Nomao.Y)]
#sample <- sample(1:24620,9844)
#Nomao.X <- Nomao.X[c(sample,24621:34464),]
#Nomao.Y <- Nomao.Y[c(sample,24621:34464)]



                                

#Ozone <- read.table("ozone.data",sep = ",")
#Ozone.X <- data.matrix(Ozone[,2:73])
#Ozone.Y <- 1+as.numeric(Ozone[,74])




#TrainX <- read.table("gisette_train.data")
#ValX <- read.table("gisette_valid.data")
#Gis.X <- rbind(as.matrix(TrainX), as.matrix(ValX))
#TrainY <- read.table("gisette_train.labels")
#ValY <- read.table("gisette_valid.labels")
#Gis.Y <- c(1/2*(3-as.numeric(TrainY[,1])),1/2*(3-as.numeric(ValY[,1])))
q

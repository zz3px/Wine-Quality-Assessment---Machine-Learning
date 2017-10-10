library(MASS)
library(e1071) #svm
wine<- read.csv("C:/Users/sony/Desktop/2016fall/STAT5330/project 2/Wine_data.csv", header = TRUE)
# scale the data
set.seed(37)
is.factor(wine$quality)
wine$quality <- as.numeric(wine$quality)

# split testing and training
t <- sample(1:nrow(wine),round(0.75*nrow(wine)))
train<-wine[t,]
test<-wine[-t,]
table(wine$quality)
dim(train)
dim(test)
# Best subset selection
library(leaps) 
regfit.best <- regsubsets(quality~.-citric.acid-chlorides-total.sulfur.dioxide,train,nvmax =11)
regfit.summary<-summary(regfit.best)


## plot
par(mfrow=c(2,2))

plot(regfit.summary$rss, xlab = "num var.",ylab = "RSS",type="l")

plot(regfit.summary$adjr2, xlab = "num var.",ylab = "adj R2",type="l")
maxpoint<-which.max(regfit.summary$adjr2) #8
points(maxpoint,regfit.summary$adjr2[maxpoint], col = "red",cex=2, pch = 20) 

plot(regfit.summary$cp, xlab = "num var.",ylab = "cp",type="l")
minpoint<-which.min(regfit.summary$cp) #8
points(minpoint,regfit.summary$cp[minpoint], col = "red",cex=2, pch = 20) 

plot(regfit.summary$bic, xlab = "num var.",ylab = "bic",type="l")
minpoint2<-which.min(regfit.summary$bic) #7
points(minpoint2,regfit.summary$bic[minpoint2], col = "red",cex=2, pch = 20) 
# indicating best model
coef(regfit.best,8)

coef(regfit.best,8)

# Validation approach
regfit.best <- regsubsets(quality~.,train,nvmax = 11)
val.errors<-rep(NA,11)  # empty set to store 8 errors
# create x matrix for prediction
test.mat<-model.matrix(quality~.,test)
for (i in 1:11){
  coefi<- coef(regfit.best,i)
  pred<-test.mat[,names(coefi)]%*%coefi     #predict y using coefficients and predictors
  val.errors[i] <- mean((test$quality-pred)^2)
}
val.errors
which.min(val.errors)  # 8 is the best
# want the final model on entire dataset
regfit.best <- regsubsets(quality~.,wine,nvmax = 11)
coef(regfit.best,8)



# Stepwise Selection
# Forward, Backward and Hybrid Stepwise
# Full model
regfit.full<-lm(quality~.,train)
# Null model
regfit.null <-lm(quality~1.,train)

# Forward Stepwise
forward<-step(regfit.null,scope = list(lower=regfit.null,upper=regfit.full), data = train,direction = "forward")
test.mat <- model.matrix(quality~.,test)
coefi <- coef(forward)
pred<-test.mat[,names(coefi)]%*%coefi
val.error<-mean((test$quality-pred)^2)
val.error #0.5973154

#backward selection
backward<-step(regfit.full,data=train,direction = "backward")
test.mat<-model.matrix(quality~.,test)
coefi2<-coef(backward)
pred2<-test.mat[,names(coefi2)]%*%coefi2     #predict y using coefficients and predictors
val.errors<- mean((test$quality-pred2)^2)  #0.5973154

# hybird slection
fit <- step(regfit.null,scope=list(lower=regfit.null,upper=regfit.full),direction = "both")
test.mat<-model.matrix(quality~.,test)
coefi<-coef(fit)
pred<-test.mat[,names(coefi)]%*%coefi     #predict y using coefficients and predictors
val.errors<- mean((test$quality-pred)^2)

#final best method on entire data
final <- lm(quality~fixed.acidity + volatile.acidity + residual.sugar 
            + free.sulfur.dioxide + density + pH + sulphates +alcohol,data = wine)


# Shrinkage
library(MASS)
library(glmnet)
#lm(y~x data)
x<-as.matrix(train[,1:11])
y<-as.matrix(train[,12])
fullx<-as.matrix(wine[,1:11])
fully<-as.matrix(wine[,12])

par(mfrow=c(1,1))
# ridge regression
fit.ridge <- glmnet(x,y,alpha = 0)
# test average error
coef<-coef(fit.ridge)
n<-ncol(coef)
aveerror<-rep(NA,n)

for (i in 1:n){
  coefi<-coef[,i]
  pred<-test.mat[,names(coefi)]%*%coefi
  aveerror[i]<-mean((test$quality-pred)^2)
}
plot(x = log(fit.ridge$lambda),y=aveerror,xlab= "log lambda",ylab="test average prediction error",main="ridge regression",type="l")
# final model on entire data ridge regression
min(aveerror) #0.5969285
bestlam<-fit.ridge$lambda[which.min(aveerror)]
final.ridge<-glmnet(fullx,fully,alpha = 0,lambda = bestlam) #0.03910185
coef(final.ridge)


# lasso regression
fit.lasso <- glmnet(x,y,alpha = 1)
# test average error
coef<-coef(fit.lasso)
n<-ncol(coef)
aveerror<-rep(NA,n)
for (i in 1:n){
  coefi<-coef[,i]
  pred<-test.mat[,names(coefi)]%*%coefi
  aveerror[i]<-mean((test$quality-pred)^2)
}
plot(x = log(fit.lasso$lambda),y=aveerror,xlab= "log lambda",ylab="test average prediction error",main="lasso regression",type="l")

# final model on entire data lasso regression
min(aveerror)  #0.5973147
bestlam<-fit.lasso$lambda[which.min(aveerror)]
final.lasso<-glmnet(fullx,fully,alpha = 1,lambda = bestlam)  # 0.001222213
coef(final.lasso)



#dimension
library(pls)

# PCR
library(pls)
pcr.fit <- pcr(quality~.,data = train, scale = T)
summary(pcr.fit)
# test average error
coef<-pcr.fit$coefficients
n<-ncol(coef)
aveerror<-rep(NA,11)

for (i in 1:11){
  coefi<-coef[,1,i]
  pred<-test.mat[,names(coefi)]%*%coefi
  aveerror[i]<-mean((test$quality-pred)^2)
}
plot(x=1:11,y=aveerror,xlab="ncomponent",ylab="test average prediction error",type="l",main = "PCR")
minpoint<-which.min(aveerror) #10
points(minpoint,aveerror[minpoint], col = "red",cex=2, pch = 20) 
min(aveerror) 
pcr.fit$loadings

# final model on entire data PCR
pcr.final <- pcr(quality~.,data=wine,ncomp = which.min(aveerror))
coef(pcr.final)


# PLS (partially linear squared regression)
pls.fit <- plsr(quality~.,data = train, scale = T)
summary(pls.fit)
pls.fit$loadings

# test average error
coef<-pls.fit$coefficients
n<-ncol(coef)
aveerror<-rep(NA,11)
for (i in 1:11){
  coefi<-coef[,1,i]
  pred<-test.mat[,names(coefi)]%*%coefi
  aveerror[i]<-mean((test$quality-pred)^2)
}
min(aveerror)
plot(x=1:11,y=aveerror,type="l",xlab="ncomponent",ylab="test average prediction error",main = "Pls")
minpoint<-which.min(aveerror) #5
points(minpoint,aveerror[minpoint], col = "red",cex=2, pch = 20) 

# best final model equation PLS
pls.final <- plsr(quality~.,data=wine,ncomp = which.min(aveerror))
coef(pls.final)
pls.final$coefficients



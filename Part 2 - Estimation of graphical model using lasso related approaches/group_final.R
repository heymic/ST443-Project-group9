#ST443 - Group project - Part II

library(MASS)
library(glmnet)
library(glasso)
library(zoo)

#Function to generate THETA

rtheta=function(p)
{
  B=matrix(0,p,p)
  #We will generate first the upper triangular matrix, then we force B to be symmetrical.
  B[upper.tri(B, diag = FALSE)]=rbinom((p^2-p)/2,1,0.1)/2
  B=B+t(B)
  #We choose diag(B) so that B is positive definite
  #If a symmetric matrix is strictly diagonal dominant then it is positive definite
  #We put on the diagonal the maximum sum by rows +1
  k=max(apply(B,1,sum))+1
  diag(B)=rep(k,p)
  #We standardize B to have unit diagonal elements    
  B=B/B[1,1]
  return(B)
}

#Function - Nodewise LASSO approach (Methods 1 and 2)
#Input - X , lambda
#Output - PRED - matrix of 0,1,2

nodewise.lasso=function(x,lambda)
{
  #par.glmnet is a 3 dimensional array - its indexes are the variable x[i], the estimated values for beta[i,j]
  #the levels if lambda
  p=ncol(x)
  par.glmnet=array(0,dim=c(p,p,length(lambda)))
  
  #We regress each X[i] on the remaining variables, we save the coefficients beta[i,j] in the 3-dimensional
  #array par.glmnet
  for(i in 1:p)
  {
    fit=glmnet(x[,-i],x[,i],lambda = lambda,intercept = FALSE)
    #We put the coefficients of fit in the 3 dimensional array par.glmnet
    for(j in 1:p)
    {
      if(j<i) par.glmnet[i,j,]=coef(fit)[j+1,]
      if(j>i) par.glmnet[i,j,]=coef(fit)[j,]
    }
  }
  
  #We estimate 1 if beta[i,j]!=0
  par.glmnet[par.glmnet!=0]=1
  
  #We arrange the 3-dimensional array into a matrix
  #For each lambda we save the estimated edges in a column of pred
  pred=matrix(0,p*(p-1)/2,length(lambda))
  
  for(i in 1:length(lambda))
  {
    mat=par.glmnet[,,(length(lambda)-i+1)]+t(par.glmnet[,,(length(lambda)-i+1)])
    pred[,i]=mat[upper.tri(mat)]
  }
  return(pred)
}

#Function - Graphical LASSO approach (Method 3)
#Input - X , lambda
#Output - PRED - matrix of 0,1,2

graphical.lasso=function(x,lambda)
{
  #Function from the library glasso
  p=ncol(x)
  var.x=var(x)
  fit.glasso=glassopath(var.x,lambda,penalize.diagonal=FALSE)
  
  #For each lambda we save the estimation for the matrix theta in the column of pred
  pred=matrix(0,p*(p-1)/2,length(lambda))
  
  for (i in 1:length(lambda))
  {
    fit=fit.glasso$wi[,,i]
    pred[,i]=as.vector(fit[upper.tri(fit)])
  }
  pred[pred!=0]=1
  return(pred)
}

#Function - ERROR calculation - TPR, FPR
#Input - Prediction of the model, true edges
#Output - TPR and FPR

error.rate=function(pred,edges,lambda)
{
  tpr=rep(0,length(lambda))
  fpr=rep(0,length(lambda))
  #Since pred and edges are vector of 0 and 1, we can use a trick to estimate the true positive rate
  #and the false positive rate. See reports for more detail.
  tpr=apply(pred*edges,2,sum)/sum(edges)
  fpr=apply(pred*(1-edges),2,sum)/sum(1-edges)
  return(cbind(tpr,fpr))
}

#Function - AUROC calculation
#We use the trapezium approximation

aur=function(tpr,fpr)
{
  #tpr=tpr[!duplicated(fpr)]
  #fpr=fpr[!duplicated(fpr)]
  o=order(fpr)
  a=sum(rollmean(tpr[o],2)*diff(fpr[o]))
  return(a)
}

#This is the MAIN function with graphical output
one.time.run=function(p,n,lambda)
{
  theta=rtheta(p)
  #sigma is the inverse of theta
  sigma=solve(theta,sparse=TRUE)
  
  #Edges is the true response variable
  edges=as.vector(theta[upper.tri(sigma)])
  edges[edges!=0]=1
  
  #We generate a n random sample from a multivariate gaussian distribution with mean 0 and variance Sigma
  x=mvrnorm(n,rep(0,p),sigma)
  
  #GLMNET package for Node-wise LASSO approach
  #Determining prediction for method 1 and 2
  #Method 1 predicts a vertex if pred=1 or pred=2
  #Method 2 predicts a vertex if pred=2
  
  pred=nodewise.lasso(x,lambda)
  #Method 1 predicts and edge if pred is equal to 1 or 2
  pred.1=matrix(0,p*(p-1)/2,length(lambda))
  pred.1[pred!=0]=1
  #Method 2 predicts and edge if pred is equal to 2
  pred.2=matrix(0,p*(p-1)/2,length(lambda))
  pred.2[pred==2]=1
  
  #GLASSO package for Graphical Lasso approach
  pred.3=graphical.lasso(x,lambda)
  
  #OUTPUT - Mis-classification rate
  #Graph
  par(mfrow=c(1,3))
  plot(lambda,2*apply(abs(pred.1-edges)/(p*(p-1)),2,sum),type="l",ylim=c(0,1),ylab = "Method 1",xlab="Lambda",col="red")
  mtext("Mis-classification rate", side = 3, line = 2, font=2)
  legend("topright",legend=c(p,n))
  plot(lambda,2*apply(abs(pred.2-edges)/(p*(p-1)),2,sum),type="l",ylim=c(0,1),ylab = "Method 2",xlab="Lambda",col="blue")
  plot(lambda,2*apply(abs(pred.3-edges)/(p*(p-1)),2,sum),type="l",ylim=c(0,1),ylab = "Method 3",xlab="Lambda",col="green")
  par(mfrow=c(1,1))
  
  #ROC and AUROC
  #positive.rate is a matrix cointaining the vector of true and false positive rate for each lambda
  
  positive.rate.1=error.rate(pred.1,edges,lambda)
  positive.rate.2=error.rate(pred.2,edges,lambda)
  positive.rate.3=error.rate(pred.3,edges,lambda)
  
  auroc.1=aur(positive.rate.1[,1],positive.rate.1[,2])
  auroc.2=aur(positive.rate.2[,1],positive.rate.2[,2])
  auroc.3=aur(positive.rate.3[,1],positive.rate.3[,2])
  
  #Plotting ROC Curves
  if(FALSE)
  {
  par(mfrow=c(1,3))
  plot(positive.rate.1[,2],positive.rate.1[,1],type="l",col="red",xlab = "False Positive Rate",ylab="True Positive Rate",ylim=c(0,1),xlim=c(0,1))
  abline(0,1,col="grey")
  text(0.6,0.2,label="AUROC = ")
  text(0.71,0.2,label=round(auroc.1,3))
  mtext("ROC curves", side = 3, line = 2, font=2)
  legend("topleft",legend=c(p,n))
  plot(positive.rate.2[,2],positive.rate.2[,1],type="l",col="blue",xlab = "False Positive Rate",ylab="True Positive Rate",ylim=c(0,1),xlim=c(0,1))
  abline(0,1,col="grey")
  text(0.6,0.2,label="AUROC = ")
  text(0.71,0.2,label=round(auroc.2,3))
  plot(positive.rate.3[,2],positive.rate.3[,1],type="l",col="green",xlab = "False Positive Rate",ylab="True Positive Rate",ylim=c(0,1),xlim=c(0,1))
  abline(0,1,col="grey")
  text(0.6,0.2,label="AUROC = ")
  text(0.71,0.2,label=round(auroc.3,3))
  par(mfrow=c(1,1))
  }
  
  li=list(list(pred.1,pred.2,pred.3),list(positive.rate.1,positive.rate.2,positive.rate.3),list(auroc.1,auroc.2,auroc.3))
  return(li) 
}

#This function is the same as before but there are not the graphical output
no.graph.run=function(p,n,lambda)
{
  theta=rtheta(p)
  sigma=solve(theta,sparse=TRUE)
  
  #Edges is the true response variable
  edges=as.vector(theta[upper.tri(sigma)])
  edges[edges!=0]=1
  
  #We generate a n random sample from a multivariate gaussian distribution with mean 0 and variance Sigma
  x=mvrnorm(n,rep(0,p),sigma)
  
  #GLMNET package for Node-wise LASSO approach
  #Determining prediction for method 1 and 2
  #Method 1 predicts a vertex if pred=1 or pred=2
  #Method 2 predicts a vertex if pred=2
  
  pred=nodewise.lasso(x,lambda)
  pred.1=matrix(0,p*(p-1)/2,length(lambda))
  pred.1[pred!=0]=1
  pred.2=matrix(0,p*(p-1)/2,length(lambda))
  pred.2[pred==2]=1
  
  #GLASSO package for Graphical Lasso approach
  pred.3=graphical.lasso(x,lambda)
  
  #ROC and AUROC
  
  positive.rate.1=error.rate(pred.1,edges,lambda)
  positive.rate.2=error.rate(pred.2,edges,lambda)
  positive.rate.3=error.rate(pred.3,edges,lambda)
  
  auroc.1=aur(positive.rate.1[,1],positive.rate.1[,2])
  auroc.2=aur(positive.rate.2[,1],positive.rate.2[,2])
  auroc.3=aur(positive.rate.3[,1],positive.rate.3[,2])

  return(cbind(auroc.1,auroc.2,auroc.3))
}

#MAIN PROGRAM
#Parameter setting

#The seed is the delivery date
set.seed(06122017)

#P=number of nodes
#n=dimension of the sample

#We define the vector lambda for the penality term
lambda=seq(0.0,0.2,by=0.002)

#result is list of list. it follow how to recall the single elements
#result.p.n[[i]][[j]]
#j=1,2,3 is the method used
#i=1 => prediction
#i=2 => true positive rate and false negative rate
#i=3 => auroc

result.50.500=one.time.run(50,500,lambda)
result.100.200=one.time.run(100,200,lambda)
result.100.500=one.time.run(100,500,lambda)
result.100.1000=one.time.run(100,1000,lambda)
result.150.500=one.time.run(150,500,lambda)

#Plot - ROC curve across methods, p=100, n=500
plot(x = c(0, 1),y = c(0, 1),type = "n",main = "ROC for method 1,2,3, p=100 , n=500",xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(0,1,col="grey")
lines(result.100.500[[2]][[1]][,2],result.100.500[[2]][[1]][,1],col="red")
lines(result.100.500[[2]][[2]][,2],result.100.500[[2]][[2]][,1],col="blue")
lines(result.100.500[[2]][[3]][,2],result.100.500[[2]][[3]][,1],col="green")

#Plot - ROC curve across methods, p=100, n=1000
plot(x = c(0, 1),y = c(0, 1),type = "n",main = "ROC for method 1,2,3, p=100 , n=500",xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(0,1,col="grey")
lines(result.100.1000[[2]][[1]][,2],result.100.1000[[2]][[1]][,1],col="red")
lines(result.100.1000[[2]][[2]][,2],result.100.1000[[2]][[2]][,1],col="blue")
lines(result.100.1000[[2]][[3]][,2],result.100.1000[[2]][[3]][,1],col="green")

#Plot - ROC curves
#Same n, p changes
par(mfrow=c(1,3))
plot(x = c(0, 1),y = c(0, 1),type = "n",main = "ROC, p=50,100,150 , n=500",xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(0,1,col="grey")
lines(result.50.500[[2]][[1]][,2],result.50.500[[2]][[1]][,1],col="red")
lines(result.100.500[[2]][[1]][,2],result.100.500[[2]][[1]][,1],col="blue")
lines(result.150.500[[2]][[1]][,2],result.150.500[[2]][[1]][,1],col="green")
plot(x = c(0, 1),y = c(0, 1),type = "n",xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(0,1,col="grey")
lines(result.50.500[[2]][[2]][,2],result.50.500[[2]][[2]][,1],col="red")
lines(result.100.500[[2]][[2]][,2],result.100.500[[2]][[2]][,1],col="blue")
lines(result.150.500[[2]][[2]][,2],result.150.500[[2]][[2]][,1],col="green")
plot(x = c(0, 1),y = c(0, 1),type = "n", xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(0,1,col="grey")
lines(result.50.500[[2]][[3]][,2],result.50.500[[2]][[3]][,1],col="red")
lines(result.100.500[[2]][[3]][,2],result.100.500[[2]][[3]][,1],col="blue")
lines(result.150.500[[2]][[3]][,2],result.150.500[[2]][[3]][,1],col="green")
par(mfrow=c(1,1))

#Plot - ROC curves
#Same p, n changes
par(mfrow=c(1,3))
plot(x = c(0, 1),y = c(0, 1),type = "n",main = "ROC, p=100 , n=200,500,1000",xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(0,1,col="grey")
lines(result.100.200[[2]][[1]][,2],result.100.200[[2]][[1]][,1],col="red")
lines(result.100.500[[2]][[1]][,2],result.100.500[[2]][[1]][,1],col="blue")
lines(result.100.1000[[2]][[1]][,2],result.100.1000[[2]][[1]][,1],col="green")
plot(x = c(0, 1),y = c(0, 1),type = "n",xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(0,1,col="grey")
lines(result.100.200[[2]][[2]][,2],result.100.200[[2]][[2]][,1],col="red")
lines(result.100.500[[2]][[2]][,2],result.100.500[[2]][[2]][,1],col="blue")
lines(result.100.1000[[2]][[2]][,2],result.100.1000[[2]][[2]][,1],col="green")
plot(x = c(0, 1),y = c(0, 1),type = "n",xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(0,1,col="grey")
lines(result.100.200[[2]][[3]][,2],result.100.200[[2]][[3]][,1],col="red")
lines(result.100.500[[2]][[3]][,2],result.100.500[[2]][[3]][,1],col="blue")
lines(result.100.1000[[2]][[3]][,2],result.100.1000[[2]][[3]][,1],col="green")
par(mfrow=c(1,1))

p=100
n=500
set.seed(06122017)
auroc=matrix(0,50,3) 
for(t in 1:50)
{
  auroc[t,]=no.graph.run(p,n,lambda)
}

apply(auroc,2,mean)
apply(auroc,2,sd)
apply(auroc,2,var)
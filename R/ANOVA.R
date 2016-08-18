#' VN ANOVA
#'
#' Analysis of variance (ANOVA) based on lower partial moment CDFs for multiple variables.  Returns a degree of certainty the samples belong to the same population, not a p-value.
#' @param A Matrix of variables.
#' @param pairwise Returns pairwise certainty tests when set to TRUE.  Defaults to FALSE.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' A<-cbind(x,y)
#' VN.ANOVA(A,pairwise=TRUE)
#' mean(na.omit(VN.ANOVA(A,pairwise = TRUE)$x))
#' @export


VN.ANOVA<- function(A,pairwise=FALSE){


  mean.of.means <- mean(colMeans(A))
  n<- ncol(A)
if(pairwise==FALSE){
  LPM_ratio = numeric(0L)
  MAD.CDF = numeric(0L)


  for (i in 1:n){
  #Continuous CDF for each variable from Mean of Means
  LPM_ratio[i] <- LPM(1,mean.of.means,A[,i])/
    (LPM(1,mean.of.means,A[,i])+UPM(1,mean.of.means,A[,i]))


  #Continuous CDF Deviation from 0.5
  MAD.CDF[i]<- abs(LPM_ratio[i] - 0.5)
  }


  Mean.MAD.CDF <- mean(MAD.CDF)


  #Certainty associated with samples
  VN.ANOVA.rho <- (0.5 - Mean.MAD.CDF)/0.5



  #Graphs
  boxplot(A, las=2, xlab= "Means", ylab="Variable",horizontal = TRUE,
           main="ANOVA", col=rainbow(n))


  #For ANOVA Visualization
  abline(v=mean.of.means,col="red",lwd=4)


  return(c("Certainty of Same Population"=VN.ANOVA.rho))

} else {

certainties = data.frame()

for(j in 1:(n-1)){
  for(i in 1:(n-j)){
    certainties[i+j,i]=VN.ANOVA(cbind(A[,i],A[,i+j]))
    certainties[i,i+j]=VN.ANOVA(cbind(A[,i],A[,i+j]))
  }
}
colnames(certainties) = colnames(A)[1:n-1]
rownames(certainties) = colnames(A)
boxplot(A, las=2, xlab= "Means", ylab="Variable",horizontal = TRUE,
        main="ANOVA", col=rainbow(n))
abline(v=mean.of.means,col="red",lwd=4)
return(certainties)
}
}

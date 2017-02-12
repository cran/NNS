#' NNS ANOVA
#'
#' Analysis of variance (ANOVA) based on lower partial moment CDFs for multiple variables.  Returns a degree of certainty the samples belong to the same population, not a p-value.
#' @param A Matrix of variables.
#' @param pairwise Returns pairwise certainty tests when set to \code{pairwise=TRUE}.  Defaults to \code{pairwise=FALSE}.
#' @return Returns the degree certainty the samples belong to the same population [0,1].
#' @keywords ANOVA
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' A<-cbind(x,y)
#' NNS.ANOVA(A,pairwise=TRUE)
#' mean(na.omit(NNS.ANOVA(A,pairwise = TRUE)$x))
#' @export


NNS.ANOVA<- function(A,pairwise=FALSE){

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
        NNS.ANOVA.rho <- (0.5 - Mean.MAD.CDF)^2/0.25

  #Graphs
        boxplot(A, las=2, xlab= "Means", ylab="Variable",horizontal = TRUE,
           main="NNS ANOVA", col=rainbow(n))

  #For ANOVA Visualization
        abline(v=mean.of.means,col="red",lwd=4)
        mtext("Grand Mean", side = 3,col = "red")

        return(c("Certainty"=NNS.ANOVA.rho))
    } else {

          certainties = data.frame()
          for(j in 1:(n-1)){
              for(i in 1:(n-j)){
                  output=NNS.ANOVA(cbind(A[,i],A[,i+j]))
                  certainties[i+j,i]=output
                  certainties[i,i+j]=output
              }
          }

          colnames(certainties) = colnames(A)
          rownames(certainties) = colnames(A)
          boxplot(A, las=2, xlab= "Means", ylab="Variable",horizontal = TRUE,
              main="ANOVA", col=rainbow(n))
          abline(v=mean.of.means,col="red",lwd=4)
          mtext("Grand Mean", side = 3,col = "red")

          return(certainties)
      }
}

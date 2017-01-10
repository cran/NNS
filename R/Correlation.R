#' NNS Correlation
#'
#' Returns the nonlinear correlation coefficient based on partial moment quadrants measured by frequency or area.  Degree = 0 is frequency, degree = 1 is area.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param order Controls the level of quadrant partitioning.  Defualts to \code{order=2}.  Errors can generally be rectified by setting \code{order=1}.
#' @param degree Defaults to 0 for smaller number of observations.
#' @return Returns nonlinear correlation coefficient between two variables, or correlation matrix for matrix input.
#' @keywords correlation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' ## Pairwise Correlation
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.cor(x,y)
#'
#' ## Correlation Matrix
#' x<-rnorm(100); y<-rnorm(100); z<-rnorm(100)
#' B<-cbind(x,y,z)
#' NNS.cor(B)
#'
#' @export

NNS.cor = function( x, y, order = 2,
                   degree= ifelse(length(x)<100,0,1)){



  if(!missing(y)){

  return(NNS.dep(x,y,print.map = F,order=order)$Correlation)

}



if(missing(y)){
  n= ncol(x)
  if(is.null(n)){stop("supply both 'x' and 'y' or a matrix-like 'x'")}
  rhos = data.frame()

  for(j in 0:(n-1)){
    for(i in 1:(n-j)){

      if((i+j)==i){
        rhos[i+j,i]=1
        rhos[i,i+j]=1
      } else {

      rhos[i+j,i]=NNS.dep(x[,i],x[,i+j],print.map = F,order=order)$Correlation
      rhos[i,i+j]=NNS.dep(x[,i],x[,i+j],print.map = F,order=order)$Correlation
      }
    }
  }
  colnames(rhos) = colnames(x)
  rownames(rhos) = colnames(x)

  return(rhos)
}


}

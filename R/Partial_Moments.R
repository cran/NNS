#' Lower Partial Moment
#'
#' This function generates a univariate lower partial moment for any degree or target.
#' @param degree \code{degree = 0} is frequency, \code{degree = 1} is area.
#' @param target Typically set to mean, but does not have to be.
#' @param variable Variable
#' @return LPM of variable
#' @keywords partial moments, mean, variance, CDF
#' @importFrom grDevices adjustcolor rainbow
#' @importFrom graphics abline boxplot legend lines par plot points segments text matplot title axis mtext
#' @importFrom stats coef cor lm na.omit sd median complete.cases resid uniroot aggregate density
#' @importFrom rgl plot3d points3d shade3d cube3d
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' LPM(0,mean(x),x)
#' @export

LPM<- function(degree,target,variable)
 {sum((target - (variable[variable <= target]))^degree)/length(variable)}



#' Upper Partial Moment
#'
#' This function generates a univariate upper partial moment for any degree or target.
#' @param degree \code{degree = 0} is frequency, \code{degree = 1} is area.
#' @param target Typically set to mean, but does not have to be.
#' @param variable Variable
#' @return UPM of variable
#' @keywords partial moments, mean, variance, upper CDF
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' UPM(0,mean(x),x)
#' @export


UPM<- function(degree,target,variable){
  sum(((variable[variable > target]) - target)^degree)/length(variable)}


#' Co-Upper Partial Moment
#' (Upper Right Quadrant 1)
#'
#' This function generates a multivariate co-upper partial moment for any degree or target.
#' @param degree.x Degree for variable X.  \code{degree.x = 0} is frequency, \code{degree.x = 1} is area.
#' @param degree.y Degree for variable Y.  \code{degree.y = 0} is frequency, \code{degree.y = 1} is area.
#' @param x Variable X
#' @param y Variable Y
#' @param target.x Typically the mean of Variable X, but does not have to be.
#' @param target.y Typically the mean of Variable Y, but does not have to be.
#' @return Co-UPM of two variables
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' Co.UPM(0,0,x,y,mean(x),mean(y))
#' @export


Co.UPM<- function(degree.x,degree.y,x,y,target.x=mean(x),target.y=mean(y)){
  x=x-target.x;y=y-target.y
  x[x<=0]<- 0;y[y<=0]<- 0
  x[x>0]<- x[x>0]^degree.x
  y[y>0]<- y[y>0]^degree.y
  return(sum(x*y)/length(x))
  }

#' Co-Lower Partial Moment
#' (Lower Left Quadrant 4)
#'
#' This function generates a multivariate co-lower partial moment for any degree or target.
#' @param degree.x Degree for variable X.  \code{degree.x = 0} is frequency, \code{degree.x = 1} is area.
#' @param degree.y Degree for variable Y.  \code{degree.y = 0} is frequency, \code{degree.y = 1} is area.
#' @param x Variable X
#' @param y Variable Y
#' @param target.x Typically the mean of Variable X, but does not have to be.
#' @param target.y Typically the mean of Variable Y, but does not have to be.
#' @return Co-LPM of two variables
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' Co.LPM(0,0,x,y,mean(x),mean(y))
#' @export

Co.LPM<- function(degree.x,degree.y,x,y,target.x=mean(x),target.y=mean(y)){
  x=target.x-x;y=target.y-y
  x[x<0]<- 0;y[y<0]<- 0
  x[x>0]<- x[x>0]^degree.x
  y[y>0]<- y[y>0]^degree.y
  return(sum(x*y)/length(x))
  }

#' Divergent-Lower Partial Moment
#' (Lower Right Quadrant 3)
#'
#' This function generates a multivariate divergent lower partial moment for any degree or target.
#' @param degree.x Degree for variable X.  \code{degree.x = 0} is frequency, \code{degree.x = 1} is area.
#' @param degree.y Degree for variable Y.  \code{degree.y = 0} is frequency, \code{degree.y = 1} is area.
#' @param x Variable X
#' @param y Variable Y
#' @param target.x Typically the mean of Variable X, but does not have to be.
#' @param target.y Typically the mean of Variable Y, but does not have to be.
#' @return Divergent LPM of two variables
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' D.LPM(0,0,x,y,mean(x),mean(y))
#' @export

D.LPM<- function(degree.x,degree.y,x,y,target.x=mean(x),target.y=mean(y)){
  x=x-target.x;y=target.y-y
  x[x<=0]<- 0;y[y<0]<- 0
  x[x>0]<- x[x>0]^degree.x
  y[y>0]<- y[y>0]^degree.y
  return(sum(x*y)/length(x))
  }


#' Divergent-Upper Partial Moment
#' (Upper Left Quadrant 2)
#'
#' This function generates a multivariate divergent upper partial moment for any degree or target.
#' @param degree.x Degree for variable X.  \code{degree.x = 0} is frequency, \code{degree.x = 1} is area.
#' @param degree.y Degree for variable Y.  \code{degree.y = 0} is frequency, \code{degree.y = 1} is area.
#' @param x Variable X
#' @param y Variable Y
#' @param target.x Typically the mean of Variable X, but does not have to be.
#' @param target.y Typically the mean of Variable Y, but does not have to be.
#' @return Divergent UPM of two variables
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' D.UPM(0,0,x,y,mean(x),mean(y))
#' @export

D.UPM<- function(degree.x,degree.y,x,y,target.x=mean(x),target.y=mean(y)){
  x=target.x-x;y=y-target.y
  x[x<0]<- 0;y[y<=0]<- 0
  x[x>0]<- x[x>0]^degree.x
  y[y>0]<- y[y>0]^degree.y
  return(sum(x*y)/length(x))
 }


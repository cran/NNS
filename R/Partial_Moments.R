#' Lower Partial Moment
#'
#' This function generates a univariate lower partial moment for any degree or target.
#' @param degree Degree = 0 is frequency, degree = 1 is area
#' @param target Typically set to mean, but does not have to be
#' @param variable Variable
#' @return LPM of variable
#' @keywords partial moments, mean, variance, CDF
#' @importFrom grDevices adjustcolor rainbow
#' @importFrom graphics abline boxplot legend lines par plot points segments text matplot title axis
#' @importFrom stats coef cor lm na.omit sd median complete.cases resid uniroot aggregate density
#' @importFrom rgl plot3d points3d
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
#' @param degree Degree = 0 is frequency, degree = 1 is area
#' @param target Typically set to mean, but does not have to be
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
#' @param degree Degree = 0 is frequency, degree = 1 is area
#' @param variable1 Variable 1
#' @param variable2 Variable 2
#' @param target1 Defaults to mean of Variable 1, but does not have to be...
#' @param target2 Defualts to mean of Variable 2, but does not have to be...
#' @return Co-UPM of two variables
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' Co.UPM(0,mean(x),mean(y),x,y)
#' @export


Co.UPM<- function(degree,variable1,variable2,target1=mean(variable1),target2=mean(variable2)){
  x=variable1-target1;y=variable2-target2
  x[x<=0]<- 0;y[y<=0]<- 0
  x[x>0]<- x[x>0]^degree
  y[y>0]<- y[y>0]^degree
  return(sum(x*y)/length(x))
  }

#' Co-Lower Partial Moment
#' (Lower Left Quadrant 4)
#'
#' This function generates a multivariate co-lower partial moment for any degree or target.
#' @param degree Degree = 0 is frequency, degree = 1 is area
#' @param variable1 Variable 1
#' @param variable2 Variable 2
#' @param target1 Defaults to mean of Variable 1, but does not have to be...
#' @param target2 Defualts to mean of Variable 2, but does not have to be...
#' @return Co-LPM of two variables
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' Co.LPM(0,mean(x),mean(y),x,y)
#' @export

Co.LPM<- function(degree,variable1,variable2,target1=mean(variable1),target2=mean(variable2)){
  x=target1-variable1;y=target2-variable2
  x[x<0]<- 0;y[y<0]<- 0
  x[x>0]<- x[x>0]^degree
  y[y>0]<- y[y>0]^degree
  return(sum(x*y)/length(x))
  }

#' Divergent-Lower Partial Moment
#' (Lower Right Quadrant 3)
#'
#' This function generates a multivariate divergent lower partial moment for any degree or target.
#' @param degree_n Degree = 0 is frequency, degree = 1 is area
#' @param degree_q Degree = 0 is frequency, degree = 1 is area
#' @param variable1 Variable 1
#' @param variable2 Variable 2
#' @param target1 Defaults to mean of Variable 1, but does not have to be...
#' @param target2 Defualts to mean of Variable 2, but does not have to be...
#' @return Divergent LPM of two variables
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' D.LPM(0,0,mean(x),mean(y),x,y)
#' @export

D.LPM<- function(degree_n,degree_q,variable1,variable2,target1=mean(variable1),target2=mean(variable2)){
  x=variable1-target1;y=target2-variable2
  x[x<=0]<- 0;y[y<0]<- 0
  x[x>0]<- x[x>0]^degree_n
  y[y>0]<- y[y>0]^degree_q
  return(sum(x*y)/length(x))
  }


#' Divergent-Upper Partial Moment
#' (Upper Left Quadrant 2)
#'
#' This function generates a multivariate divergent upper partial moment for any degree or target.
#' @param degree_n Degree = 0 is frequency, degree = 1 is area
#' @param degree_q Degree = 0 is frequency, degree = 1 is area
#' @param variable1 Variable 1
#' @param variable2 Variable 2
#' @param target1 Defaults to mean of Variable 1, but does not have to be...
#' @param target2 Defualts to mean of Variable 2, but does not have to be...
#' @return Divergent UPM of two variables
#' @keywords partial moments, covariance
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' D.UPM(0,0,mean(x),mean(y),x,y)
#' @export

D.UPM<- function(degree_n,degree_q,variable1,variable2,target1=mean(variable1),target2=mean(variable2)){
  x=target1-variable1;y=variable2-target2
  x[x<0]<- 0;y[y<=0]<- 0
  x[x>0]<- x[x>0]^degree_n
  y[y>0]<- y[y>0]^degree_q
  return(sum(x*y)/length(x))
 }


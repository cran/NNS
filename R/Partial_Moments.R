#' Lower Partial Moment
#'
#' This function generates a univariate lower partial moment for any degree or target.
#' @param degree Degree = 0 is frequency, degree = 1 is area
#' @param target Typically set to mean, but does not have to be
#' @param variable Variable
#' @return LPM of variable
#' @keywords partial moments
#' @importFrom grDevices adjustcolor rainbow
#' @importFrom graphics abline boxplot legend lines par plot points segments text
#' @importFrom stats coef cor lm na.omit sd
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' \dontrun{LPM(0,0,x)}
#' @export

LPM<- function(degree,target,variable)
 {sum((target - (variable[variable < target]))^degree)/length(variable)}



#' Upper Partial Moment
#'
#' This function generates a univariate upper partial moment for any degree or target.
#' @param degree Degree = 0 is frequency, degree = 1 is area
#' @param target Typically set to mean, but does not have to be
#' @param variable Variable
#' @return UPM of variable
#' @keywords partial moments
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' \dontrun{UPM(0,0,x)}
#' @export


UPM<- function(degree,target,variable){
  sum(((variable[variable > target]) - target)^degree)/length(variable)}


#' Co-Upper Partial Moment
#' (Upper Right Quadrant 1)
#'
#' This function generates a multivariate co-upper partial moment for any degree or target.
#' @param degree Degree = 0 is frequency, degree = 1 is area
#' @param target1 Typically set to mean of Variable 1, but does not have to be
#' @param target2 Typically set to mean of Variable 2, but does not have to be
#' @param variable1 Variable 1
#' @param variable2 Variable 2
#' @return Co-UPM of two variables
#' @keywords partial moments
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{Co_UPM(0,mean(x),mean(y),x,y)}
#' @export


Co_UPM<- function(degree,target1,target2,variable1,variable2){
  output <- vector("numeric", length(variable1))
  for (i in 1:length(variable1))
  {
    if ((variable1[i]>target1)*(variable2[i]>target2)==1)

      output[i]<- (((variable1[i]-target1)^degree)*((variable2[i]-target2)^degree))
  }

  return(Co_UPM=sum(output)/length(variable1))

}

#' Co-Lower Partial Moment
#' (Lower Left Quadrant 4)
#'
#' This function generates a multivariate co-lower partial moment for any degree or target.
#' @param degree Degree = 0 is frequency, degree = 1 is area
#' @param target1 Typically set to mean of Variable 1, but does not have to be
#' @param target2 Typically set to mean of Variable 2, but does not have to be
#' @param variable1 Variable 1
#' @param variable2 Variable 2
#' @return Co-LPM of two variables
#' @keywords partial moments
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{Co_LPM(0,mean(x),mean(y),x,y)}
#' @export

Co_LPM<- function(degree,target1,target2,variable1,variable2){

  output <- vector("numeric", length(variable1))
  for (i in 1:length(variable1))
  {
    if (variable1[i]<target1 & variable2[i]<target2)

    output[i]<- (((target1-variable1[i])^degree)*((target2-variable2[i])^degree))
  }

  return(sum(output)/length(variable1))
}

#' Divergent-Lower Partial Moment
#' (Lower Right Quadrant 3)
#'
#' This function generates a multivariate divergent lower partial moment for any degree or target.
#' @param degree_n Degree = 0 is frequency, degree = 1 is area
#' @param degree_q Degree = 0 is frequency, degree = 1 is area
#' @param target1 Typically set to mean of Variable 1, but does not have to be
#' @param target2 Typically set to mean of Variable 2, but does not have to be
#' @param variable1 Variable 1
#' @param variable2 Variable 2
#' @return Divergent LPM of two variables
#' @keywords partial moments
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{D_LPM(0,mean(x),mean(y),x,y)}
#' @export

D_LPM<- function(degree_n,degree_q,target1,target2,variable1,variable2){

  output <- vector("numeric", length(variable1))
  for (i in 1:length(variable1))
  {


    if ((variable1[i]>target1)*(variable2[i]<target2)==1)

      output[i]<- (((variable1[i]-target1)^degree_q)*((target2-variable2[i])^degree_n))
  }
  return(c(sum(output)/length(variable1)))
}


#' Divergent-Upper Partial Moment
#' (Upper Left Quadrant 2)
#'
#' This function generates a multivariate divergent upper partial moment for any degree or target.
#' @param degree_n Degree = 0 is frequency, degree = 1 is area
#' @param degree_q Degree = 0 is frequency, degree = 1 is area
#' @param target1 Typically set to mean of Variable 1, but does not have to be
#' @param target2 Typically set to mean of Variable 2, but does not have to be
#' @param variable1 Variable 1
#' @param variable2 Variable 2
#' @return Divergent UPM of two variable
#' @keywords partial moments
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{D_UPM(0,mean(x),mean(y),x,y)}
#' @export

D_UPM<- function(degree_n,degree_q,target1,target2,variable1,variable2){

  output <- vector("numeric", length(variable1))
  for (i in 1:length(variable1))
  {
    if ((variable1[i]<target1)*(variable2[i]>target2)==1)

      output[i]<- (((target1-variable1[i])^degree_q)*((variable2[i]-target2)^degree_n))
  }
  return(c(sum(output)/length(variable1)))
}


#' Nonlinear FWL
#'
#' Applies \code{VN.reg} to residuals per FWL (Frisch-Waugh-Lovell) Theorem for multiple nonlinear regression coefficients.
#'
#'
#' @param B Complete dataset of independent variables (IV) in matrix form.
#' @param y Dependent variable (DV).
#' @param order Controls the number of the \code{VN.reg}.  Defaults to \code{order=1}.
#' @param type Controls the partitioning in \code{VN.reg}.  Set to \code{type="XONLY"} for IV based partitioning.  Defaults to NULL for both IV and DV partitioning.
#' @param linear.test Set to \code{linear.test=TRUE} to run linear regression on residuals and compare to multiple linear regression output.
#' @keywords Frisch-Waugh-Lovell, multiple nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123);size=30;x_1=rnorm(size,mean=5,sd=1);x_2=rnorm(size,mean=1,sd=5); x_3=runif(size)
#' y=x_1^2+x_2+x_3
#' B=cbind(x_1,x_2,x_3)
#' NNS.FWL(B,y)
#'
#' ## Test the coefficients
#' set.seed(123);size=30;x_1=rnorm(size,mean=5,sd=1);x_2=rnorm(size,mean=1,sd=5); x_3=runif(size)
#' y=x_1^2+x_2+x_3
#' B=cbind(x_1,x_2,x_3)
#' NNS.FWL(B,y,linear.test=TRUE)
#' @export


NNS.FWL <- function (B,y,order=1,type=NULL,linear.test=FALSE){
  n=ncol(B)


  ### For single regressor, Detrend X and Y by regressing each on the time variable.
  if(is.null(n)){
    time.index = 1:length(B)
    x.star=numeric()
    y.star=numeric()
    x.star = mean(B)+resid(lm(B~time.index))
    y.star = mean(y)+resid(lm(y~time.index))

    if(linear.test==FALSE){
      print("Raw Value VN.reg")
      print(VN.reg(x,y,return.equation = TRUE,type = type,order=order)$derivative)
      print("Linear Residuals VN.reg")
      return((VN.reg(x.star,y.star,order=order,return.equation = TRUE,type = type)$derivative))
    }

    if(linear.test==TRUE){
      print("Linear Residuals Linear Regression")
      print(coef(lm(y~B)))
      return(coef(lm(y.star~x.star)))
    }

  }

  ### For Multiple regressions
  if(!is.null(n)){
    original.variable = B
    new.variable = matrix(nrow=nrow(B))

    ###  Turn each column into numeric values
    for (i in 1:ncol(original.variable)){

      new.variable = cbind(new.variable,as.numeric(original.variable[,i]))
    }

    B=new.variable[,-1]
    colnames(B)=colnames(original.variable)
    y=as.numeric(y)

    #  Regress x^k_i on all K-1 variables and store residuals
    r=list()
    ry=list()
    Beta=list()

    #  Name and define all residual vectors r and ry
    for(i in 1:n){
      r[[paste('r',i,sep="")]] <- numeric()
      r[[i]]<- resid(lm(B[,i]~B[,-i]))
      ry[[paste('ry',i,sep="")]] <- numeric()
      ry[[i]]<- resid(lm(y~B[,-i]))

      if(linear.test==FALSE){
        ## FOR NNS Beta
        derivs <- VN.reg(r[[i]],ry[[i]],return.equation = TRUE,order = order,plot = FALSE,type = type)$derivative[,1]
        lower.x <- VN.reg(B[,i],y,return.equation = TRUE,order = order,plot = FALSE,type = type)$derivative[,2]
        upper.x <- VN.reg(B[,i],y,return.equation = TRUE,order = order,plot = FALSE,type = type)$derivative[,3]

        Beta[[paste('Beta_',i,sep="")]] <- data.frame(matrix(nrow=length(derivs)))
        Beta[[i]] <- matrix(c(derivs,lower.x,upper.x),byrow=FALSE,nrow=length(derivs))
        colnames(Beta[[i]])<- c('Coefficient','X Lower Range','X Upper Range')
      }

      if(linear.test==TRUE){
        ## LINEAR BETA TEST
        Beta[[paste('Beta_',i,sep="")]] <- data.frame(matrix(ncol=n))
        Beta[[i]] <- coef(lm(ry[[i]]~r[[i]]))[2]
      }
    }

  }



if(linear.test==FALSE){
  p = length((Beta[[1]][,1]))


  y.hat.matrix = matrix(ncol=n,nrow=length(y))
  y.hat.list=list()

  for(j in 1:n){
  x=B[,j]
  fitted = matrix(ncol=2)
  fitted.new = numeric()

    for (i in 1:p){

        z=(which(x>=Beta[[j]][i,2] & x<=Beta[[j]][i,3]))

        z.diff = ((x[z]- Beta[[j]][i,2])*Beta[[j]][i,1])+Beta[[j]][i,2]

        fitted.new =  cbind(z,z.diff)


        fitted= rbind(fitted,fitted.new)
        fitted = fitted[order(fitted[,1]),]
       }


  y.hat.list[[j]]=fitted[,2]

    }

    y.hat.matrix=do.call('cbind',y.hat.list)

    y.hat=(na.omit(rowSums(y.hat.matrix)))


    R2=  (sum((y.hat-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((y.hat-mean(y))^2))
}

    print("Raw Values Multiple Linear Regression")
    print(coef(lm(y~B)))
    print("Regression on Residuals")
    VN.reg.output = VN.reg(B,y,order=order,type = type,return.values = TRUE)$derivative
    return(c(Beta,VN.reg.output))

}

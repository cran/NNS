#' VN Classifier
#'
#' Classifies data based on multiple logistic nonlinear regression \code{VN.reg}
#'
#' @param x Complete dataset in matrix form.
#' @param y Column of data to be classified.
#' @param threshold  Sets the correlation threshold for independent variables.  Defaults to 0.
#' @param order Controls the number of partial moment quadrant means.  Users are encouraged to try different \code{order=} settings.  \code{order='max'} will use maximum suggested possible order based on number of observations.
#' @param point.est Returns the fitted value for any value of the independent variable.  Use a vector of values for independent varaiables to return the multiple regression fitted value.
#' @param location Sets the legend location within the plot
#' @param print.values Defaults to FALSE, set to TRUE in order to return all fitted values for independent variable
#' @param print.equation Defaults to FALSE, set to TRUE in order to return the local coefficients
#' @param plot  To plot regression or not.  Defaults to TRUE.
#' @param clean.method  Method to handle missing or NA values.  'omit' uses \code{\link{complete.cases}} function on data while 'zero' replaces missing data with 0 value.  Defaults to NULL assuming user has cleaned data prior to analyzing.
#' @param dep.order Sets the internal order for \link{VN.dep}.  Categorical variables typically require \code{dep.order=1}.  Error message will alert user if this is the case.
#' @return Returns two variables, mean squared error "\code{MSE}" and fitted values "\code{Fitted}" as well as the number of predictors, R2, R2 Adjusted, and Prediction Accuracy measured by percentage of exact classifications.
#' @keywords nonlinear logistic regression, classifier
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' ## Using 'iris' dataset where predictive attributes are columns 1:4, and the class is column 5.
#' VN.class(iris,5)
#'
#' ## To call mean squared error
#' VN.class(iris,5)$MSE
#'
#' ## To call fitted values
#' VN.class(iris,5)$Fitted
#' @export



VN.class = function (x, y,threshold = 0,
                     order=NULL,
                     point.est = NULL,
                     location = 'top',
                     print.values = FALSE,
                     print.equation = FALSE,
                     plot = TRUE,
                     clean.method = NULL,
                     dep.order=NULL){


  type = "CLASS"

  preds = numeric()
  x[is.na(x)] <- 0
  original.columns = ncol(x)
  original.variable = x
  new.variable = matrix(nrow=nrow(x))


  ###  Turn each column into numeric values
  for (i in 1:ncol(original.variable)){

    new.variable = cbind(new.variable,as.numeric(original.variable[,i]))

    ###  Clean data
    if(!is.null(clean.method)){
      if(clean.method == 'omit'){
        new.variable = new.variable[complete.cases(new.variable)]
      }

      if(clean.method == 'zero'){
        new.variable[is.na(new.variable)] <- 0
      }
    }
  }

  x <- new.variable[, c(-1, -(y + 1))]
  y <- new.variable[, (y + 1)]
  total.cor = Co.PM.cor(cbind(x,y))


  if(!is.null(ncol(x))) {

    if (ncol(x)==1){
      x=(x)
    }else{
      x.star.dep = numeric()
      x.star.coef = numeric()
      x.star.matrix = matrix(nrow=length(y))

      for (i in 1:ncol(x)){

        x.star.dep[i] = VN.dep(x[,i],y,print.map=FALSE,order=dep.order)[1]

        x.star.coef[i] = x.star.dep[i]-total.cor
        if(is.na(x.star.coef[i])){return("CATEGORICAL VARIABLE, PLEASE SET dep.order=1 ")}
        if(abs(x.star.coef[i])<threshold){x.star.coef[i]=0}

        x.star.matrix =  cbind(x.star.matrix,x.star.coef[i]*x[,i])


        if(i == ncol(x)){

          print(paste0("Synthetic Independent Variable X* = (",
                       paste(format(x.star.coef[1:i],digits = 4),paste("X",1:i,sep = ''),sep='*',collapse = "  "),")/",sum(abs(x.star.coef)>0)),quote=FALSE)

        }
      }


      if(!is.null(point.est)){

        point.est= sum(point.est*x.star.coef)/sum(abs(x.star.coef)>0)
      } else {point.est=NULL}


      x = rowSums(x.star.matrix[,2:(1+ncol(x))])/ncol(x)
    }}



  if(is.null(order)){

    reg.output = VN.reg(x,y,type = "CLASS",return.values = TRUE,plot = FALSE,point.est = point.est)
    regression.points = reg.output$regression.points
    Regression.Coefficients= reg.output$derivative

    y.fitted = reg.output$Fitted[,2]


    if(!is.null(point.est)){
      point.est.y=reg.output$Point.est
    } else{point.est.y=NULL}

    if(print.equation==TRUE){
      print(Regression.Coefficients)
    }
    Fitted.values = y.fitted

    Values = (cbind(x,Fitted=y.fitted,Actual=y,Difference=y.fitted-(y),
                    Accuracy=abs(round(y.fitted)-(y))
    ))

    MSE = mean((y.fitted-y)^2)

    R2=(sum((y.fitted-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((y.fitted-mean(y))^2))

    R2.adj = R2

    Prediction.Accuracy=(length(y)-sum(abs(round(y.fitted)-(y))>0))/length(y)

    preds=  MSE
  }




  if(!is.null(order)){
    if(order=="max"){
      order=ceiling(log2(length(y)))
    }else{order = order}

    reg.output = VN.reg(x,y,type = "CLASS",order = order,return.values = TRUE,plot = FALSE,point.est=point.est)

    regression.points = reg.output$regression.points
    Regression.Coefficients= reg.output$derivative
    y.fitted = reg.output$Fitted[,2]
    if(!is.null(point.est)){
      point.est.y=reg.output$Point.est
    } else{point.est.y=NULL}

    if(print.equation==TRUE){
      print(Regression.Coefficients)
    }

    Values = (cbind(x,Fitted=y.fitted,Actual=y,Difference=y.fitted-(y),
                    Accuracy=abs(round(y.fitted)-(y))
    ))



    MSE = mean((y.fitted-y)^2)

    Fitted.values = y.fitted

    R2=(sum((y.fitted-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((y.fitted-mean(y))^2))

    R2.adj = R2

    Prediction.Accuracy=(length(y)-sum(abs(round(y.fitted)-(y))>0))/length(y)

  }

  p = length(reg.output$derivative[,2])

  ###Plotting and regression equation
  if(plot==TRUE){
    xmin= min(c(point.est,x))
    xmax= max(c(point.est,x))
    ymin= min(c(point.est.y,y))
    ymax= max(c(point.est.y,y))
    plot(x,y,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),col='steelblue',
         xlab = if(original.columns>1){"Synthetic X*"}else{"X"},
         ylab="Y",main=paste0("Order = ",order))


    ### Plot Regression points and fitted values and legend
    points(na.omit(regression.points[order(regression.points),]),col='red',pch=19)
    lines(na.omit(regression.points[order(regression.points),]),col='red',lwd=2,lty = 2)


    if(!is.null(point.est)){ points(point.est,point.est.y, col='green',pch=18)
      legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),
                                                          paste("Segments",(p-1)),
                                                          if(original.columns>1){paste("Synthetic Point Estimate",point.est)}else{paste("Point Estimate",point.est)},
                                                          paste("Fitted Value",format(point.est.y,digits = 6))
      ))}

    if(is.null(point.est)){
      legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),paste("Segments",(p-1))))
    }

    if(!is.null(point.est)){
      if(point.est>max(x)) segments(point.est,point.est.y,regression.points[p+ceiling(abs(p-q)/2),1],regression.points[p+ceiling(abs(p-q)/2),2],col="green",lty=2)
      if(point.est<min(x)) segments(point.est,point.est.y,regression.points[1,1],regression.points[1,2],col="green",lty=2)
    }

  }# plot TRUE bracket



  ### Print Values

  if(print.values ==FALSE){
    if(is.null(point.est)){
      print(c("Segments" = p-1,"R2"=R2,
              Prediction.Accuracy=Prediction.Accuracy

      ))
      return(list("MSE"=MSE,"Fitted"=as.numeric(Fitted.values)))
    }}

  if(print.values ==FALSE){
    if(!is.null(point.est)) {
      print(c("Segments" = p-1,"R2"=R2))
      if(!is.null(original.columns)){if(original.columns>1){
        print(c(Synthetic_Point=point.est, Fitted.value=point.est.y
        ))
        return(list("MSE"=MSE,"Fitted"=as.numeric(Fitted.values)))}
      }else{ print(c(Point=point.est, Fitted.value=point.est.y))
        return(list("MSE"=MSE,"Fitted"=as.numeric(Fitted.values)))}
    }}

  if(print.values ==TRUE){
    if(is.null(point.est)){
      print(Values)
      print(c("Segments" = p-1,"R2"=R2,
              Prediction.Accuracy=Prediction.Accuracy
      ))
      return(list("MSE"=MSE,"Fitted"=as.numeric(Fitted.values)))
    }}


  if(print.values ==TRUE){
    if(!is.null(point.est)) {
      print(Values)
      print(c("Segments" = p-1,"R2"=R2))
      if(!is.null(original.columns)){if(original.columns>1){
        print(c(Synthetic_Point=point.est, Fitted.value=point.est.y))
        return(list("MSE"=MSE,"Fitted"=as.numeric(Fitted.values)))}
      }else{ print(c(Point=point.est, Fitted.value=point.est.y))
        return(list("MSE"=MSE,"Fitted"=as.numeric(Fitted.values)))}
    }}



}

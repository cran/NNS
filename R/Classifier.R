#' VN Classifier
#'
#' Classifies data based on multiple logistic nonlinear regression (VN.reg)
#'
#' @param x Complete dataset in matrix form.
#' @param y Column of data to be classified.
#' @param threshold  Sets the coefficient threshold for independent variables.  Defaults to 0.
#' @param order Controls the number of partial moment quadrant means.  Defaults to smaller order to avoid overfitting
#' @param point.est Returns the fitted value for any value of the independent variable.  Use a vector of values for independent varaiables to return the multiple regression fitted value.
#' @param location Sets the legend location within the plot
#' @param print.values Defaults to FALSE, set to TRUE in order to return all fitted values for independent variable
#' @param print.equation Defaults to FALSE, set to TRUE in order to return the local coefficients
#' @param plot  To plot regression or not.  Defaults to TRUE.
#' @keywords nonlinear logistic regression, classifier
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' ## Using 'iris' dataset where predictive attributes are columns 1:4, and the class is column 5.
#' \dontrun{VN.class(iris,5)}
#' @export


VN.class = function (x, y,threshold = 0,
                     order=NULL,
                     point.est = NULL,
                     location = 'top',
                     print.values = FALSE,
                     print.equation = FALSE,
                     plot = TRUE){

  type = "LOGIT"

  preds = numeric()

  original.columns = ncol(x)
  original.variable = x
  new.variable = matrix(nrow=nrow(x))

  ###  Turn each column into numeric values
  for (i in 1:ncol(original.variable)){

    new.variable = cbind(new.variable,as.numeric(original.variable[,i]))

    ###  Clean data
    new.variable = new.variable[,!colSums(!is.finite(new.variable))]
  }

  x= new.variable[,-y]
  y= new.variable[,y]


  if(!is.null(ncol(x))) {

    if (ncol(x)==1){
      x=(x)
    }else{

      x.star.coef = numeric()
      x.star.matrix = matrix(nrow=length(y))



        for (i in 1:ncol(x)){


          for (j in 2:floor(log(length(x),4))){
            if(is.na(VN.cor(x[,i],y,j,degree=0))){
              cor.order=j-1
              break}
          }



        x.star.coef[i] = VN.cor(x[,i],y,cor.order)

        if(abs(x.star.coef[i])<threshold){x.star.coef[i]=0}

        x.star.matrix =  cbind(x.star.matrix,x.star.coef[i]*x[,i])



        if(i == ncol(x)){

          print(paste0("Synthetic Independent Variable X* = (",
                       paste(format(x.star.coef[1:i],digits = 4),paste("X",1:i,sep = ''),sep='*',collapse = "  "),")/",sum(abs(x.star.coef)>0)),quote=FALSE)
          print("",quote=FALSE)
        }
      }


      if(!is.null(point.est)){

        point.est= sum(point.est*x.star.coef)/sum(abs(x.star.coef)>0)
      }


      x = rowSums(x.star.matrix[,2:(1+ncol(x))])/ncol(x)
    }}

if(is.null(order)){

for (k in 1:ceiling(log(length(y),2))){

   order.p=k

  temp_df = data.frame(x=x, y=y)
  temp_df[,'temp_part'] = 'p'
  temp_df[,'master_part'] = 'p'

  regression.points = data.frame(matrix(ncol = 2))
  Regression.Coefficients = data.frame(matrix(ncol=3))

  names(Regression.Coefficients) = c('Coefficient','X Lower Range','X Upper Range')



 # if(order<1){return("Please Increase the Order Specification")}

  ###  Partition Map based on X-values only
  if(!is.null(type)){
    order.p=order.p-1

    for(i in 0:(order.p)){

      for(item in unique(temp_df$master_part)){
        tmp_xbar = mean(temp_df[temp_df$master_part == item,'x'])
        tmp_ybar = mean(temp_df[temp_df$master_part == item,'y'])



        temp_df[temp_df$x >= tmp_xbar  & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar  & temp_df$master_part == item,'master_part'], 1, sep = '')
        temp_df[temp_df$x < tmp_xbar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x < tmp_xbar  & temp_df$master_part == item,'master_part'], 2, sep = '')



        ### order + 1 to account for 'p'
        if(nchar(item)==order.p+1){

          regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

        }


      }

      temp_df[,'master_part'] = temp_df[, 'temp_part']

    }


  }



  x0 = unique(y[x==min(x)])
  x.max = unique(y[x==max(x)])

  if(length(x0)>1){x0 = mean(x0)}
  if(length(x.max)>1){x.max = mean(x.max)}


  regression.points[1,2] = x0
  regression.points[1,1] = min(x)

  regression.points[length(regression.points[,2])+1,2] = x.max
  regression.points[length(regression.points[,1]),1] = max(x)


  ###Regression Equation

  regression.points = na.omit(regression.points[order(regression.points),])


  q=length(regression.points[,1])



  for(i in 1:q){

    rise = regression.points[i+1,2] - regression.points[i,2]
    run = regression.points[i+1,1] - regression.points[i,1]

    Regression.Coefficients[i,] = cbind((rise/run),regression.points[i,1],regression.points[i+1,1])
    Regression.Coefficients[q,] = cbind(1,regression.points[i,1],regression.points[i,1]+1e-10)
  }

  Regression.Coefficients= na.omit(Regression.Coefficients)

  ### Fitted Values
  p = length((Regression.Coefficients)[,1])

  ### Differences in p, q arise at times...
  if(p!=q){break}

  fitted = numeric()
  fitted.new = numeric()

  for (i in 1:p){

    z=(which(x>=Regression.Coefficients[i,2] & x<Regression.Coefficients[(i),3]))

    z.diff = ((x[z]- Regression.Coefficients[i,2])*Regression.Coefficients[i,1])+regression.points[i,2]


    if(is.null(point.est)){point.est.y = NULL} else{

      if(!is.null(point.est) && point.est>=Regression.Coefficients[i,2] && point.est<Regression.Coefficients[i,3]){ point.est.y = (point.est - Regression.Coefficients[i,2])*(Regression.Coefficients[i,1])+regression.points[i,2]}

      else{if(!is.null(point.est) && point.est<Regression.Coefficients[1,2]){
        point.est.y = ((point.est - Regression.Coefficients[1,2])*(Regression.Coefficients[1,1]))+(regression.points[1,2])
      }

        else{if(!is.null(point.est) && point.est>Regression.Coefficients[p,2]){point.est.y = ((point.est - Regression.Coefficients[(p-0),2])*(Regression.Coefficients[(p-1),1]))+(regression.points[(p),2])
        }
        }
      }
    }


    fitted.new =  cbind(z,z.diff)


    fitted = rbind(fitted,fitted.new)
    fitted = fitted[order(fitted[,1]),]

  }

  if(print.equation==TRUE){
    print(regression.points)
    print(Regression.Coefficients)
  }

  Values = (cbind(x,Fitted=fitted[,2],Actual=y,Difference=fitted[,2]-(y),
                  Accuracy=abs(round(fitted[,2])-(y))
  ))

  MSE = mean((fitted[,2]-y)^2)



  R=cor(fitted[,2],y)
  R2=R^2

  R2.adj = 1 - (((1-R2)*length(fitted))/(length(fitted)-p-1))

  Prediction.Accuracy=(length(y)-sum(abs(round(fitted[,2])-(y))>0))/length(y)

  ###Plotting and regression equation
  if(plot==TRUE){
    xmin= min(c(point.est,x))
    xmax= max(c(point.est,x))
    ymin= min(c(point.est.y,y))
    ymax= max(c(point.est.y,y))
    plot(x,y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),col='steelblue',
         xlab = if(original.columns>1){"Synthetic X*"}else{"X"},
         ylab="Y",main=paste0("Order = ",order+1))


    ### Plot Regression points and fitted values and legend
    points(na.omit(regression.points[order(regression.points),]),col='red',pch=19)
    lines(na.omit(regression.points[order(regression.points),]),col='red',lwd=2,lty = 2)


    if(!is.null(point.est)){ points(point.est,point.est.y, col='green',pch=18)
      legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),
                                                          paste("R2 Adjusted",format(R2.adj,digits=4)),paste("Predictors",(p-1)),
                                                          if(original.columns>1){paste("Synthetic Point Estimate",point.est)}else{paste("Point Estimate",point.est)},
                                                          paste("Fitted Value",format(point.est.y,digits = 6))
      ))}

    if(is.null(point.est)){
      legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),paste("R2 Adjusted",format(R2.adj,digits=4)),paste("Predictors",(p-1))))
    }

    if(!is.null(point.est)){
      if(point.est>max(x)) segments(point.est,point.est.y,regression.points[p,1],regression.points[p,2],col="green",lty=2)
      if(point.est<min(x)) segments(point.est,point.est.y,regression.points[1,1],regression.points[1,2],col="green",lty=2)
    }

  }# plot TRUE bracket
  preds[k]=  Prediction.Accuracy

} # k order bracket
}

  optimal.order = which.max(na.omit(preds))


  if(!is.null(order)){order = order} else {order=optimal.order}
  #order=optimal.order

  temp_df = data.frame(x=x, y=y)
  temp_df[,'temp_part'] = 'p'
  temp_df[,'master_part'] = 'p'

  regression.points = data.frame(matrix(ncol = 2))
  Regression.Coefficients = data.frame(matrix(ncol=3))

  names(Regression.Coefficients) = c('Coefficient','X Lower Range','X Upper Range')



  # if(order<1){return("Please Increase the Order Specification")}

  ###  Partition Map based on X-values only
  if(!is.null(type)){
    order=order-1

    for(i in 0:(order)){

      for(item in unique(temp_df$master_part)){
        tmp_xbar = mean(temp_df[temp_df$master_part == item,'x'])
        tmp_ybar = mean(temp_df[temp_df$master_part == item,'y'])



        temp_df[temp_df$x >= tmp_xbar  & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar  & temp_df$master_part == item,'master_part'], 1, sep = '')
        temp_df[temp_df$x < tmp_xbar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x < tmp_xbar  & temp_df$master_part == item,'master_part'], 2, sep = '')



        ### order + 1 to account for 'p'
        if(nchar(item)==order+1){

          regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

        }


      }

      temp_df[,'master_part'] = temp_df[, 'temp_part']

    }


  }



  x0 = unique(y[x==min(x)])
  x.max = unique(y[x==max(x)])

  if(length(x0)>1){x0 = mean(x0)}
  if(length(x.max)>1){x.max = mean(x.max)}


  regression.points[1,2] = x0
  regression.points[1,1] = min(x)

  regression.points[length(regression.points[,2])+1,2] = x.max
  regression.points[length(regression.points[,1]),1] = max(x)


  ###Regression Equation

  regression.points = na.omit(regression.points[order(regression.points),])


  q=length(regression.points[,1])



  for(i in 1:q){

    rise = regression.points[i+1,2] - regression.points[i,2]
    run = regression.points[i+1,1] - regression.points[i,1]

    Regression.Coefficients[i,] = cbind((rise/run),regression.points[i,1],regression.points[i+1,1])
    Regression.Coefficients[q,] = cbind(1,regression.points[i,1],regression.points[i,1]+1e-10)
  }

  Regression.Coefficients= na.omit(Regression.Coefficients)

  ### Fitted Values
  p = length((Regression.Coefficients)[,1])

  ### Differences in p, q arise at times...
  if(p!=q){break}

  fitted = numeric()
  fitted.new = numeric()

  for (i in 1:p){

    z=(which(x>=Regression.Coefficients[i,2] & x<Regression.Coefficients[(i),3]))

    z.diff = ((x[z]- Regression.Coefficients[i,2])*Regression.Coefficients[i,1])+regression.points[i,2]


    if(is.null(point.est)){point.est.y = NULL} else{

      if(!is.null(point.est) && point.est>=Regression.Coefficients[i,2] && point.est<Regression.Coefficients[i,3]){ point.est.y = (point.est - Regression.Coefficients[i,2])*(Regression.Coefficients[i,1])+regression.points[i+ceiling(abs(p-q)/2),2]}

      else{if(!is.null(point.est) && point.est<Regression.Coefficients[1,2]){
        point.est.y = ((point.est - Regression.Coefficients[1,2])*(Regression.Coefficients[1,1]))+(regression.points[1,2])
      }

        else{if(!is.null(point.est) && point.est>Regression.Coefficients[p,2]){point.est.y = ((point.est - Regression.Coefficients[(p-0),2])*(Regression.Coefficients[(p-1),1]))+(regression.points[(p+ceiling(abs(p-q)/2)),2])
        }
        }
      }
    }


    fitted.new =  cbind(z,z.diff)


    fitted = rbind(fitted,fitted.new)
    fitted = fitted[order(fitted[,1]),]

  }

  if(print.equation==TRUE){
    #print(regression.points)
    print(Regression.Coefficients[-p,])
  }

  Values = (cbind(x,Fitted=fitted[,2],Actual=y,Difference=fitted[,2]-(y),
                  Accuracy=abs(round(fitted[,2])-(y))
  ))

  MSE = mean((fitted[,2]-y)^2)


  R2=(sum((fitted[,2]-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((fitted[,2]-mean(y))^2))

  R2.adj = 1 - (((1-R2)*length(fitted))/(length(fitted)-p-1))

  Prediction.Accuracy=(length(y)-sum(abs(round(fitted[,2])-(y))>0))/length(y)

  ###Plotting and regression equation
  if(plot==TRUE){
    xmin= min(c(point.est,x))
    xmax= max(c(point.est,x))
    ymin= min(c(point.est.y,y))
    ymax= max(c(point.est.y,y))
    plot(x,y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),col='steelblue',
         xlab = if(original.columns>1){"Synthetic X*"}else{"X"},
         ylab="Y",main=paste0("Order = ",order+1))


    ### Plot Regression points and fitted values and legend
    points(na.omit(regression.points[order(regression.points),]),col='red',pch=19)
    lines(na.omit(regression.points[order(regression.points),]),col='red',lwd=2,lty = 2)


    if(!is.null(point.est)){ points(point.est,point.est.y, col='green',pch=18)
      legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),
                                                          paste("R2 Adjusted",format(R2.adj,digits=4)),paste("Predictors",(p-1)),
                                                          if(original.columns>1){paste("Synthetic Point Estimate",point.est)}else{paste("Point Estimate",point.est)},
                                                          paste("Fitted Value",format(point.est.y,digits = 6))
      ))}

    if(is.null(point.est)){
      legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),paste("R2 Adjusted",format(R2.adj,digits=4)),paste("Predictors",(p-1))))
    }

    if(!is.null(point.est)){
      if(point.est>max(x)) segments(point.est,point.est.y,regression.points[p+ceiling(abs(p-q)/2),1],regression.points[p+ceiling(abs(p-q)/2),2],col="green",lty=2)
      if(point.est<min(x)) segments(point.est,point.est.y,regression.points[1,1],regression.points[1,2],col="green",lty=2)
    }

  }# plot TRUE bracket


























































    ### Print Values

  if(print.values ==FALSE){
    if(is.null(point.est)){
      return(c("Predictors" = p-1,"R2"=R2,"R2 Adjusted"=R2.adj,
               Prediction.Accuracy=(length(y)-sum(abs(round(fitted[,2])-(y))>0))/length(y)

      ))
    }}

  if(print.values ==FALSE){
    if(!is.null(point.est)) {
      print(c("Predictors" = p-1,"R2"=R2,"R2 Adjusted"=R2.adj))
      if(!is.null(original.columns)){if(original.columns>1){
        return(c(Synthetic_Point=point.est, Fitted.value=point.est.y
        ))}
      }else{ return(c(Point=point.est, Fitted.value=point.est.y))}
    }}

  if(print.values ==TRUE){
    if(is.null(point.est)){
      print(Values)
      return(c("Predictors" = p-1,"R2"=R2,"R2 Adjusted"=R2.adj,
               Prediction.Accuracy=(length(y)-sum(abs(round(fitted[,2])-(y))>0))/length(y)
      ))
    }}


  if(print.values ==TRUE){
    if(!is.null(point.est)) {
      print(Values)
      print(c("Predictors" = p-1,"R2"=R2,"R2 Adjusted"=R2.adj))
      if(!is.null(original.columns)){if(original.columns>1){
        return(c(Synthetic_Point=point.est, Fitted.value=point.est.y))}
      }else{ return(c(Point=point.est, Fitted.value=point.est.y))}
    }}



}

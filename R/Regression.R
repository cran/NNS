#' VN Regression
#'
#' Generates a nonlinear regression based on partial moment quadrant means.
#'
#' @param x Independent Variable(s)
#' @param y Dependent Variable
#' @param order Controls the number of partial moment quadrant means.  Defaults to smaller order to avoid overfitting
#' @param type  To perform logistic regression, set to type = "LOGIT".  Defualts to NULL.
#' @param point.est Returns the fitted value for any value of the independent variable.  Use a vector of values for independent varaiables to return the multiple regression fitted value.
#' @param location Sets the legend location within the plot
#' @param print.values Defaults to FALSE, set to TRUE in order to print all fitted values for independent variable
#' @param print.equation Defaults to FALSE, set to TRUE in order to print the local coefficients (partial derivative wrt independen variable) for a given range of the independent variable
#' @param return.values Defaults to FALSE, set to TRUE in order to return fitted values into a vector.
#' @param return.equation Defaults to FALSE, set to TRUE in order to return the local coefficients (partial derivative wrt independen variable) for a given range of the independent variable
#' @param plot  To plot regression or not.  Defaults to TRUE.
#' @param tol  Minimum acceptable R2 Adjusted to stop order testing.  Defaults to 0.99999.
#' @keywords nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.reg(x,y)}
#'
#' ## For Multiple Regression:
#' x<-cbind(rnorm(100),rnorm(100),rnorm(100)); y<-rnorm(100)
#' \dontrun{VN.reg(x,y,point.est=c(.25,.5,.75))}
#'
#' ## To call fitted values:
#' \dontrun{VN.reg(x,y,return.values=TRUE)$fitted}
#'
#' ## To call partial derivative:
#' \dontrun{VN.reg(x,y,return.equation=TRUE)$derivative}
#'
#' @export


VN.reg = function (x, y,
                   order=NULL,
                   type = NULL,
                   point.est = NULL,
                   location = 'top',
                   print.values = FALSE,
                   print.equation = FALSE,
                   return.values = FALSE,
                   return.equation = FALSE,
                   plot = TRUE,
                   tol=.99999){

  R2s = numeric()
  original.columns = ncol(x)
  original.variable = x

if(is.null(order)){
for (k in 1:ceiling(log(length(y),2))){

  order.p=k



  if(!is.null(ncol(original.variable))) {
if (ncol(original.variable)==1){
      x=(original.variable)
    }else{


      x.star.coef = numeric()
      x.star.matrix = matrix(nrow=length(y))



      for (i in 1:ncol(original.variable)){


    #    for (j in 1:floor(log(length(y),2))){
     #     if(is.na(VN.dep(original.variable[,i],y,j,print.map=FALSE)[1])){
     #       cor.order=j-1
     #       break}
     #   }

        cor.order=1

        x.star.coef[i]=  VN.dep(original.variable[,i],y,cor.order,print.map = FALSE)[1]
        x.star.matrix =  cbind(x.star.matrix,x.star.coef[i]*original.variable[,i])


      }


      if(!is.null(point.est)){

        point.est= sum(point.est*x.star.coef)/ncol(original.variable)
      }


      x = rowSums(x.star.matrix[,2:(1+ncol(original.variable))])/ncol(original.variable)
    }}



  temp_df = data.frame(x=x, y=y)
  temp_df[,'temp_part'] = 'p'
  temp_df[,'master_part'] = 'p'


  regression.points = data.frame(matrix(ncol = 2))
  Regression.Coefficients = data.frame(matrix(ncol=3))

  colnames(Regression.Coefficients) = c('Coefficient','X Lower Range','X Upper Range')

 # if(order<1){return("Please Increase the Order Specification")}
  if(is.null(type)){
order.p = order.p-1
    for(i in 0:(order.p)){

      for(item in unique(temp_df$master_part)){


        tmp_xbar = mean(temp_df[temp_df$master_part == item,'x'])
        tmp_ybar = mean(temp_df[temp_df$master_part == item, 'y'])



        temp_df[temp_df$x >= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'master_part'], 1, sep = '')
        temp_df[temp_df$x < tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x < tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'master_part'], 2, sep = '')
        temp_df[temp_df$x >= tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'master_part'], 3, sep = '')
        temp_df[temp_df$x < tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x < tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'master_part'], 4, sep = '')


        ### order + 1 to account for 'p'
        if(nchar(item)==order.p+1){

          regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

        }


      }

      temp_df[,'master_part'] = temp_df[, 'temp_part']

    }
  }


  ###  FOR LOGISTIC REGRESSION
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


  for (j in 1:ceiling(log(length(y),2))){
    if(is.na(VN.dep(x,y,j,print.map=FALSE)[2])){
      dep.order=j-1
      break}
  }



  min.range = min(na.omit(regression.points[,1]))
  max.range = max(na.omit(regression.points[,1]))


  Dynamic.average.min = median(y[x<min.range])
  Dynamic.average.max = median(y[x>max.range])

  ###Endpoints
  if(is.null(type)){
    if(length(x[x<min.range])>0){
      if(VN.dep(x,y,dep.order,print.map = FALSE)[2]<.5){
        x0 = Dynamic.average.min} else {
          x0 = unique(y[x==min(x)])} }  else {x0 = unique(y[x==min(x)])}

    if(length(x[x>max.range])>0){
      if(VN.dep(x,y,dep.order,print.map = FALSE)[2]<.5){x.max = Dynamic.average.max} else {x.max = unique(y[x==max(x)])}}  else { x.max = unique(y[x==max(x)])}
  }

  if(!is.null(type)){
    x0 = unique(y[x==min(x)])
    x.max = unique(y[x==max(x)])

    if(length(x0)>1){x0 = mean(x0)}
    if(length(x.max)>1){x.max = mean(x.max)}

  }

  regression.points[1,2] = x0
  regression.points[1,1] = min(x)

  regression.points[length(regression.points[,2])+1,2] = x.max
  regression.points[length(regression.points[,1]),1] = max(x)



  ###Regression Equation

  regression.points = na.omit(regression.points[order(regression.points),])


  q=length(regression.points[,2])



  for(i in 1:q){

    rise = regression.points[i+1,2] - regression.points[i,2]
    run = regression.points[i+1,1] - regression.points[i,1]

    Regression.Coefficients[i,] = cbind((rise/run),regression.points[i,1],regression.points[i+1,1])
    Regression.Coefficients[q,] = cbind(1,regression.points[i,1],regression.points[i,1]+1e-10)
  }

  Regression.Coefficients= na.omit(Regression.Coefficients)

  ### Fitted Values
  p = length((Regression.Coefficients)[,2])

  ### Differences in p, q arise at times...
  if(p!=q){break}

  fitted = numeric()
  fitted.new = numeric()


  for (i in 1:p){

    z=(which(x>=Regression.Coefficients[i,2] & x<Regression.Coefficients[i,3]))

    z.diff = ((x[z]- Regression.Coefficients[i,2])*Regression.Coefficients[i,1])+regression.points[i,2]


    if(is.null(point.est)){point.est.y = NULL} else{

      if(point.est>=Regression.Coefficients[i,2] & point.est<Regression.Coefficients[i,3]){

        point.est.y = (point.est - Regression.Coefficients[i,2])*(Regression.Coefficients[i,1])+regression.points[i,2]
        }

      else{if(point.est<Regression.Coefficients[1,2]){
        point.est.y = ((point.est - Regression.Coefficients[1,2])*(Regression.Coefficients[1,1]))+(regression.points[1,2])
      }

        else{if(point.est>Regression.Coefficients[p,2]){point.est.y = ((point.est - Regression.Coefficients[(p-0),2])*(Regression.Coefficients[(p-1),1]))+(regression.points[(p),2])
        }
        }
      }
    }


    fitted.new =  cbind(z,z.diff)


    fitted = rbind(fitted,fitted.new)
    fitted = fitted[order(fitted[,1]),]

  }



  Values = (cbind(x,Fitted=fitted[,2],Actual=y,Difference=fitted[,2]-(y),
                  Accuracy=abs(round(fitted[,2])-(y))
  ))

  MSE = mean((fitted[,2]-y)^2)



  R=cor(fitted[,2],y)

  R2=  (sum((fitted[,2]-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((fitted[,2]-mean(y))^2))

  R2.adj = 1 - (((1-R2)*length(fitted))/(length(fitted)-p-1))

  ###Plotting and regression equation
  if(plot==TRUE){
    xmin= min(c(point.est,x))
    xmax= max(c(point.est,x))
    ymin= min(c(point.est.y,y))
    ymax= max(c(point.est.y,y))
    plot(x,y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),col='steelblue',
         xlab = if(!is.null(original.columns))
         {if(original.columns>1){"Synthetic X*"}}else{"X"},
         ylab="Y",main=paste0("Order = ",order+1))



    ### Plot Regression points and fitted values and legend
    points(na.omit(regression.points[order(regression.points),]),col='red',pch=19)
    lines(na.omit(regression.points[order(regression.points),]),col='red',lwd=2,lty = 2)


    if(!is.null(point.est)){ points(point.est,point.est.y, col='green',pch=18)
      legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),
                                                          paste("R2 Adjusted",format(R2.adj,digits=4)),paste("Predictors",(p-1)),
                                                          if(!is.null(original.columns))
                                                          {if(original.columns>1){paste("Synthetic Point Estimate",point.est)}}else{paste("Point Estimate",point.est)},
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

R2s[k]=  R2.adj

if(!is.null(tol)){
if(R2s[k]>tol) break}
}}
optimal.order = which.max(na.omit(R2s))


###  FOR RUNNING WITH OPTIMAL ORDER
if(!is.null(order)){order = order} else {order=optimal.order}

if(!is.null(ncol(original.variable))) {

  if (ncol(original.variable)==1){
    x=(original.variable)
  }else{

    x.star.coef = numeric()
    x.star.matrix = matrix(nrow=length(y))



    for (i in 1:ncol(original.variable)){



    #  for (j in 1:floor(log(length(y),2))){
    #    if(is.na(VN.cor(original.variable[,i],y,j))){
    #      cor.order=j-1
    #      break}
    #  }

      cor.order = 1

      x.star.coef[i]=  VN.dep(original.variable[,i],y,cor.order,print.map = FALSE)[1]
      x.star.matrix =  cbind(x.star.matrix,x.star.coef[i]*original.variable[,i])



      if(i == ncol(original.variable)){

        print(paste0("Synthetic Independent Variable X* = (",
                     paste(format(x.star.coef[1:i],digits = 4),paste("X",1:i,sep = ''),sep='*',collapse = "  "),")/",ncol(original.variable)),quote=FALSE)
        print("",quote=FALSE)
      }
    }


    if(!is.null(point.est)){

      point.est= sum(point.est*x.star.coef)/ncol(original.variable)
    }


    x = rowSums(x.star.matrix[,2:(1+ncol(original.variable))])/ncol(original.variable)
  }}



temp_df = data.frame(x=x, y=y)
temp_df[,'temp_part'] = 'p'
temp_df[,'master_part'] = 'p'

regression.points = data.frame(matrix(ncol = 2))
Regression.Coefficients = data.frame(matrix(ncol=3))

colnames(Regression.Coefficients) = c('Coefficient','X Lower Range','X Upper Range')


if(is.null(type)){
order=order-1

  for(i in 0:(order)){

    for(item in unique(temp_df$master_part)){


      tmp_xbar = mean(temp_df[temp_df$master_part == item,'x'])
      tmp_ybar = mean(temp_df[temp_df$master_part == item, 'y'])



      temp_df[temp_df$x >= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'master_part'], 1, sep = '')
      temp_df[temp_df$x < tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x < tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'master_part'], 2, sep = '')
      temp_df[temp_df$x >= tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'master_part'], 3, sep = '')
      temp_df[temp_df$x < tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x < tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'master_part'], 4, sep = '')


      ### order + 1 to account for 'p'
      if(nchar(item)==order+1){

        regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

      }


    }

    temp_df[,'master_part'] = temp_df[, 'temp_part']

  }
}


###  FOR LOGISTIC REGRESSION
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

for (j in 1:floor(log(length(x),2))){
  if(is.na(VN.dep(x,y,j,print.map=FALSE)[2])){
    dep.order=j-1
    break}
}

min.range = min(na.omit(regression.points[,1]))
max.range = max(na.omit(regression.points[,1]))


Dynamic.average.min = median(y[x<min.range])
Dynamic.average.max = median(y[x>max.range])

###Endpoints
if(is.null(type)){
  if(length(x[x<min.range])>0){
    if(VN.dep(x,y,dep.order,print.map = FALSE)[2]<.5){
      x0 = Dynamic.average.min} else {
        x0 = unique(y[x==min(x)])} }  else {x0 = unique(y[x==min(x)])}

  if(length(x[x>max.range])>0){
    if(VN.dep(x,y,dep.order,print.map = FALSE)[2]<.5){x.max = Dynamic.average.max} else {x.max = unique(y[x==max(x)])}}  else { x.max = unique(y[x==max(x)])}
}

if(!is.null(type)){
  x0 = unique(y[x==min(x)])
  x.max = unique(y[x==max(x)])

  if(length(x0)>1){x0 = mean(x0)}
  if(length(x.max)>1){x.max = mean(x.max)}

}

regression.points[1,2] = x0
regression.points[1,1] = min(x)

regression.points[length(regression.points[,2])+1,2] = x.max
regression.points[length(regression.points[,1]),1] = max(x)



###Regression Equation

regression.points = na.omit(regression.points[order(regression.points),])


q=length(regression.points[,2])



for(i in 1:q){

  rise = regression.points[i+1,2] - regression.points[i,2]
  run = regression.points[i+1,1] - regression.points[i,1]

  Regression.Coefficients[i,] = cbind((rise/run),regression.points[i,1],regression.points[i+1,1])
  Regression.Coefficients[q,] = cbind(1,regression.points[i,1],regression.points[i,1]+1e-10)
}

Regression.Coefficients= na.omit(Regression.Coefficients)

### Fitted Values
p = length((Regression.Coefficients)[,2])


fitted = numeric()
fitted.new = numeric()

for (i in 1:p){

  z=(which(x>=Regression.Coefficients[i,2] & x<Regression.Coefficients[(i),3]))

  z.diff = ((x[z]- Regression.Coefficients[i,2])*Regression.Coefficients[i,1])+regression.points[i,2]


  if(is.null(point.est)){point.est.y = NULL} else{

    if( point.est>=Regression.Coefficients[i,2] && point.est<Regression.Coefficients[i,3]){ point.est.y = (point.est - Regression.Coefficients[i,2])*(Regression.Coefficients[i,1])+regression.points[i,2]}

    else{if( point.est<Regression.Coefficients[1,2]){
      point.est.y = ((point.est - Regression.Coefficients[1,2])*(Regression.Coefficients[1,1]))+(regression.points[1,2])
    }

      else{if( point.est>Regression.Coefficients[p,2]){point.est.y = ((point.est - Regression.Coefficients[(p-0),2])*(Regression.Coefficients[(p-1),1]))+(regression.points[(p),2])
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


R=cor(fitted[,2],y)^2


R2=  (sum((fitted[,2]-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((fitted[,2]-mean(y))^2))

R2.adj = 1 - (((1-R2)*length(fitted))/(length(fitted)-p-1))

###Plotting and regression equation
if(plot==TRUE){
  xmin= min(c(point.est,x))
  xmax= max(c(point.est,x))
  ymin= min(c(point.est.y,y))
  ymax= max(c(point.est.y,y))
  plot(x,y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),col='steelblue',
       xlab = if(!is.null(original.columns))
       {if(original.columns>1){"Synthetic X*"}}else{"X"},
       ylab="Y",main=paste0("Order = ",order+1))



  ### Plot Regression points and fitted values and legend
  points(na.omit(regression.points[order(regression.points),]),col='red',pch=19)
  lines(na.omit(regression.points[order(regression.points),]),col='red',lwd=2,lty = 2)


  if(!is.null(point.est)){ points(point.est,point.est.y, col='green',pch=18)
    legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),
                                                        paste("R2 Adjusted",format(R2.adj,digits=4)),paste("Predictors",(p-1)),
                                                        if(!is.null(original.columns))
                                                        {if(original.columns>1){paste("Synthetic Point Estimate",point.est)}}else{paste("Point Estimate",point.est)},
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


### Print / return Values

  if(return.values == TRUE | return.equation==TRUE){
    return(list("fitted"=cbind(x,"fitted values"=fitted[,2]), "derivative"=Regression.Coefficients[-p,]))
  }

  if(print.values ==FALSE){
    if(is.null(point.est)){
      return(c("Predictors" = p-1,"R2"=R2,"R2 Adjusted"=R2.adj))
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
      return(c("Predictors" = p-1,"R2"=R2,"R2 Adjusted"=R2.adj))
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




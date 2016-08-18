#' VN Regression
#'
#' Generates a nonlinear regression based on partial moment quadrant means.
#'
#' @param x Independent Variable(s)
#' @param y Dependent Variable
#' @param order Controls the number of partial moment quadrant means.  Users are encouraged to try different \code{order=} settings.  \code{order='max'} will use maximum suggested possible order based on number of observations.
#' @param type  To perform logistic regression, set to \code{type = "LOGIT"}.  Defualts to NULL.
#' @param point.est Returns the fitted value for any value of the independent variable.  Use a vector of values for independent varaiables to return the multiple regression fitted value.
#' @param location Sets the legend location within the plot
#' @param print.values Defaults to FALSE, set to TRUE in order to print all fitted values for independent variable
#' @param print.equation Defaults to FALSE, set to TRUE in order to print the local coefficients (partial derivative wrt independen variable) for a given range of the independent variable
#' @param return.values Defaults to FALSE, set to TRUE in order to return fitted values into a vector.
#' @param return.equation Defaults to FALSE, set to TRUE in order to return the local coefficients (partial derivative wrt independent variable) for a given range of the independent variable
#' @param plot  To plot regression or not.  Defaults to TRUE.
#' @param threshold  Sets the correlation threshold for independent variables.  Defaults to 0.
#' @keywords nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' VN.reg(x,y)
#'
#' ## Manual {order} selection
#' VN.reg(x,y,order=2)
#'
#' ## Maximum {order} selection
#' VN.reg(x,y,order='max')
#'
#' ## x-only paritioning
#' VN.reg(x,y,type="XONLY")
#'
#' ## For Multiple Regression based on Synthetic X*:
#' x<-cbind(rnorm(100),rnorm(100),rnorm(100)); y<-rnorm(100)
#' VN.reg(x,y,point.est=c(.25,.5,.75))
#'
#' ## To call fitted values:
#' VN.reg(x,y,return.values=TRUE)$fitted
#'
#' ## To call partial derivative:
#' VN.reg(x,y,return.equation=TRUE)$derivative
#'
#' @export


VN.reg = function (x,y,
                   order=NULL,
                   type = NULL,
                   point.est = NULL,
                   location = 'top',
                   print.values = FALSE,
                   print.equation = FALSE,
                   return.values = FALSE,
                   return.equation = FALSE,
                   plot = TRUE,
                   threshold = 0){

  R2s = numeric()
  original.columns = ncol(x)
  original.variable = x
  total.cor = Co.PM.cor(cbind(x,y))



if(!is.null(ncol(original.variable))) {
        if (ncol(original.variable)==1){
          x=(original.variable)
        }else{

          type = "XONLY"
          x.star.dep = numeric()
          x.star.coef = numeric()
          x.star.matrix = matrix(nrow=length(y))


          for (i in 1:ncol(original.variable)){
            x.star.structure = VN.dep(original.variable[,i],y,print.map = FALSE)

            x.star.dep[i] = x.star.structure[2]
            x.star.coef[i]=  x.star.structure[1]-total.cor
            if(abs(x.star.coef[i])<threshold){x.star.coef[i]=0}
            x.star.matrix =  cbind(x.star.matrix,x.star.coef[i]*original.variable[,i])

            if(i == ncol(x)){

              print(paste0("Synthetic Independent Variable X* = (",
                           paste(format(x.star.coef[1:i],digits = 4),paste("X",1:i,sep = ''),sep='*',collapse = "  "),")/",sum(abs(x.star.coef)>0)),quote=FALSE)
              print("",quote=FALSE)
            }


          }


          if(!is.null(point.est)){

            point.est= sum(point.est*x.star.coef)/ncol(original.variable)
          }


          x = rowSums(x.star.matrix[,2:(1+ncol(original.variable))])/ncol(original.variable)
        }}



  if (is.null(original.columns)){
  dependence = VN.dep(x,y,print.map = FALSE)[2]}else{dependence=mean(x.star.dep)}

  if(is.null(order)){
      if(!is.null(original.columns)|!is.null(type)){
        if(dependence>0.75&length(y)>=100){
        part.map = partition.map(x,y,overide=TRUE)
          if(!is.null(type)){part.map=partition.map(x,y,type="XONLY",overide=TRUE)}}

        if(dependence<0.75|length(y)<100){
          part.map = partition.map(x,y)
          if(!is.null(type)){part.map=partition.map(x,y,type="XONLY")}}
      } else {

        if(dependence>0.75&length(y)>=100){
        part.map = partition.map(x,y,overide=TRUE)
          if(!is.null(type)){part.map=partition.map(x,y,overide=TRUE)}}

        if(dependence<0.75|length(y)<100){
        part.map = partition.map(x,y)
          if(!is.null(type)){part.map=partition.map(x,y)}}
      }

      Regression.Coefficients = data.frame(matrix(ncol=3))

      colnames(Regression.Coefficients) = c('Coefficient','X Lower Range','X Upper Range')

      regression.points=part.map$regression.points
      naive.order =min(nchar(part.map$df[,3]))

      min.range = min(regression.points[,1])
      max.range = max(regression.points[,1])


      Dynamic.average.min = median(y[x<min.range])
      Dynamic.average.max = median(y[x>max.range])

      ###Endpoints
      if(is.null(type)){
        if(length(x[x<min.range])>0){
          if(dependence<.5){
            x0 = Dynamic.average.min} else {
              x0 = unique(y[x==min(x)])} }  else {x0 = unique(y[x==min(x)])}

        if(length(x[x>max.range])>0){
          if(dependence<.5){x.max = Dynamic.average.max} else {x.max = unique(y[x==max(x)])}}  else { x.max = unique(y[x==max(x)])}
      }

      if(!is.null(type)){
        x0 = unique(y[x==min(x)])
        x.max = unique(y[x==max(x)])

        if(length(x0)>1){x0 = mean(x0)}
        if(length(x.max)>1){x.max = mean(x.max)}

      }

      regression.points = rbind(regression.points,c(min(x),x0))
      regression.points = rbind(regression.points,c(max(x),x.max))


### Regression Equation

      regression.points = na.omit(regression.points[order(regression.points),])


      q=length(regression.points[,2])
  ### Eliminate possible differences in p, q
      if(regression.points[q,1]==regression.points[(q-1),1]){regression.points=regression.points[-q,]}
      if(regression.points[1,1]==regression.points[2,1]){regression.points=regression.points[-1,]}

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


      R2=  (sum((fitted[,2]-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((fitted[,2]-mean(y))^2))

      R2.adj = R2 #1 - (((1-R2)*length(fitted))/(length(fitted)-p-1))

  ###Plotting and regression equation
      if(plot==TRUE){
        xmin= min(c(point.est,x))
        xmax= max(c(point.est,x))
        ymin= min(c(point.est.y,y))
        ymax= max(c(point.est.y,y))
        plot(x,y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),col='steelblue',
             xlab = if(!is.null(original.columns))
             {if(original.columns>1){"Synthetic X*"}}else{"X"},
             ylab="Y",main=paste0("Order = ",naive.order))



    ### Plot Regression points and fitted values and legend
        points(na.omit(regression.points[order(regression.points),]),col='red',pch=19)
        lines(na.omit(regression.points[order(regression.points),]),col='red',lwd=2,lty = 2)


        if(!is.null(point.est)){ points(point.est,point.est.y, col='green',pch=18)
          legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),
                                                              paste("Segments",(p-1)),
                                                              if(!is.null(original.columns))
                                                              {if(original.columns>1){paste("Synthetic Point Estimate",point.est)}}else{paste("Point Estimate",point.est)},
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
      optimal.order = min(nchar(part.map$df[,3]))-1

    }




  ###  FOR RUNNING WITH OPTIMAL ORDER
    if(!is.null(order)){
      if(order=="max"){
        order=ceiling(log2(length(y)))
      }else{order = order}

      if(is.null(type)){part.map = partition.map(x,y,overide = TRUE,order = order)}else{part.map=partition.map(x,y,type="XONLY",overide = TRUE,order = order)}


      regression.points = data.frame(matrix(ncol = 2))
      Regression.Coefficients = data.frame(matrix(ncol=3))

      colnames(Regression.Coefficients) = c('Coefficient','X Lower Range','X Upper Range')


      regression.points=part.map$regression.points

      min.range = min(regression.points[,1])
      max.range = max(regression.points[,1])

      Dynamic.average.min = median(y[x<min.range])
      Dynamic.average.max = median(y[x>max.range])

  ###Endpoints
      if(is.null(type)){
          if(length(x[x<min.range])>0){
              if(dependence<.5){
                  x0 = Dynamic.average.min} else {
                  x0 = unique(y[x==min(x)])} }  else {x0 = unique(y[x==min(x)])}

          if(length(x[x>max.range])>0){
              if(dependence<.5){x.max = Dynamic.average.max} else {x.max = unique(y[x==max(x)])}}  else { x.max = unique(y[x==max(x)])}
  }

      if(!is.null(type)){
          x0 = unique(y[x==min(x)])
          x.max = unique(y[x==max(x)])

          if(length(x0)>1){x0 = mean(x0)}
          if(length(x.max)>1){x.max = mean(x.max)}

  }

  regression.points = rbind(regression.points,c(min(x),x0))
  regression.points = rbind(regression.points,c(max(x),x.max))



  ###Regression Equation

    regression.points = na.omit(regression.points[order(regression.points),])


    q=length(regression.points[,2])

    if(regression.points[q,1]==regression.points[q-1,1]){regression.points=regression.points[-q,]}
    if(regression.points[2,1]==regression.points[1,1]){regression.points=regression.points[-1,]}

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
    print(Regression.Coefficients[-p,])
  }

  Values = (cbind(x,Fitted=fitted[,2],Actual=y,Difference=fitted[,2]-(y),
                  Accuracy=abs(round(fitted[,2])-(y))
  ))

  MSE = mean((fitted[,2]-y)^2)

  R2=  (sum((fitted[,2]-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((fitted[,2]-mean(y))^2))

  R2.adj = R2#1 - (((1-R2)*length(fitted))/(length(fitted)-p-1))

  ###Plotting and regression equation
    if(plot==TRUE){
      xmin= min(c(point.est,x))
      xmax= max(c(point.est,x))
      ymin= min(c(point.est.y,y))
      ymax= max(c(point.est.y,y))
      plot(x,y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),col='steelblue',
         xlab = if(!is.null(original.columns))
         {if(original.columns>1){"Synthetic X*"}}else{"X"},
         ylab="Y",main=paste0("Order = ",order))



    ### Plot Regression points and fitted values and legend
      points(na.omit(regression.points[order(regression.points),]),col='red',pch=19)
      lines(na.omit(regression.points[order(regression.points),]),col='red',lwd=2,lty = 2)


      if(!is.null(point.est)){ points(point.est,point.est.y, col='green',pch=18)
      legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),
                                                          paste("Segments",(p-1)),
                                                          if(!is.null(original.columns))
                                                          {if(original.columns>1){paste("Synthetic Point Estimate",point.est)}}else{paste("Point Estimate",point.est)},
                                                          paste("Fitted Value",format(point.est.y,digits = 6))
      ))}

      if(is.null(point.est)){
      legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),paste("Segments",(p-1))))
    }

      if(!is.null(point.est)){
          if(point.est>max(x)) segments(point.est,point.est.y,regression.points[p,1],regression.points[p,2],col="green",lty=2)
          if(point.est<min(x)) segments(point.est,point.est.y,regression.points[1,1],regression.points[1,2],col="green",lty=2)
    }

  }# plot TRUE bracket
}

  ### Print / return Values
    if(return.values == TRUE | return.equation==TRUE ){
        return(list("fitted"=cbind(x,"fitted values"=fitted[,2]), "derivative"=Regression.Coefficients[-p,],"Point.est"=point.est.y,"regression.points"=regression.points,"R2"=R2))
  }



    if(print.values ==FALSE){
        if(is.null(point.est)){
          return(c("Segments" = p-1,"R2"=R2))
    }}

    if(print.values ==FALSE){
        if(!is.null(point.est)) {
          print(c("Segments" = p-1,"R2"=R2))
        if(!is.null(original.columns)){if(original.columns>1){
            return(c(Synthetic_Point=point.est, Fitted.value=point.est.y
        ))}}else{ return(c(Point=point.est, Fitted.value=point.est.y))}
    }}

  if(print.values ==TRUE){
    if(is.null(point.est)){
      print(Values)
      return(c("Segments" = p-1,"R2"=R2))
    }}


  if(print.values ==TRUE){
    if(!is.null(point.est)) {
      print(Values)
      print(c("Segments" = p-1,"R2"=R2))
      if(!is.null(original.columns)){if(original.columns>1){
        return(c(Synthetic_Point=point.est, Fitted.value=point.est.y))}
      }else{ return(c(Point=point.est, Fitted.value=point.est.y))}
    }}



}

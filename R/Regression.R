#' NNS Regression
#'
#' Generates a nonlinear regression based on partial moment quadrant means.
#'
#' @param x Independent Variable(s)
#' @param y Dependent Variable
#' @param order Controls the number of partial moment quadrant means.  Users are encouraged to try different \code{order=} integer settings with \code{noise.reduction='off'}.  \code{order="max"} will force limit perfect fit.
#' @param s.t.n Signal to noise parameter, sets the threshold of \code{NNS.dep} which reduces \code{"order"} when \code{order=NULL}.  Defaults to 0.99 to ensure high dependence for higher \code{"order"} and endpoint determination.
#' @param type  To perform logistic regression, set to \code{type = "LOGIT"}.  To perform a classification, set to \code{type = "CLASS"}.  Defualts to NULL.
#' @param point.est Returns the fitted value for any value of the independent variable.  Use a vector of values for independent varaiables to return the multiple regression fitted value.
#' @param location Sets the legend location within the plot
#' @param return.values Defaults to TRUE, set to FALSE in order to only display a regression plot.
#' @param plot  To plot regression or not.  Defaults to TRUE.
#' @param residual.plot To plot the \code{fitted values of Y} and \code{Y}.  Defaults to TRUE.
#' @param threshold  Sets the correlation threshold for independent variables.  Defaults to 0.
#' @param dep.order Sets the internal order for \link{NNS.dep}.  Categorical variables typically require \code{dep.order=1}.  Error message will alert user if this is the case.  Defaults to \code{dep.order=NULL}.
#' @param n.best Sets the number of nearest regression points to use in kernel weighting for multivariate regression.  Defaults to 1.
#' @param precision  Increases speed of computation at the expense of precision.  3 settings offered: \code{"LOW"} (Default setting), \code{"MED"}, and \code{"HIGH"}.  \code{"HIGH"} is the limit condition of every observation as a regression point.
#' @param text If performing a text classification, set \code{text=TRUE}.  Defaults to FALSE.
#' @param noise.reduction In low signal:noise situations,\code{noise.reduction="mean"} (Default setting) uses means for \link{NNS.dep} restricted partitions, \code{noise.reduction="median"} uses medians instead of means for \link{NNS.dep} restricted partitions, while \code{noise.reduction="mode"}  uses modes instead of means for \link{NNS.dep} restricted partitions.  \code{noise.reduction='off'}  allows for maximum possible fit and specific \code{order} specification.
#' @param norm Normalizes regressors between 0 and 1 for multivariate regression when set to \code{norm="std"}, or normalizes regressors according to \link{NNS.norm} when set to \code{norm="NNS"}. Defaults to NULL.
#' @param dist Selects the distance calculation used. \code{dist="L2"} (default) selects the Euclidean distance and \code{dist="L1"} seclects the Manhattan distance.
#' @return UNIVARIATE regression returns the values:  \code{"Fitted"} for only the fitted values of the DV; \code{"Fitted.xy"} for a data frame of IV and fitted values; \code{"derivative"} for the coefficient of the IV and its applicable range; \code{"Point"} returns the IV point(s) being evaluated; \code{"Point.est"} for the predicted value generated; \code{"regression.points"} provides the points used in the regression equation for the given order of partitions; \code{"R2"} provides the goodness of fit.
#'
#' MULTIVARIATE regression returns the values: \code{"Fitted"} for only the fitted values of the DV; \code{"Fitted.xy"} for a data frame of IV and fitted values; \code{"regression.points"} provides the points for each IV used in the regression equation for the given order of partitions; \code{"rhs.partitions"} returns the partition points for each IV; \code{"partition"} returns the DV, quadrant assigned to the observation and fitted value; \code{"Point"} returns the IV point(s) being evaluated; \code{"Point.est"} returns the predicted value generated; \code{"equation"} returns the synthetic X* dimension reduction equation.
#' @keywords nonlinear regression, classifier
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.reg(x,y)
#'
#' ## Manual {order} selection
#' NNS.reg(x,y,order=2)
#'
#' ## Maximum {order} selection
#' NNS.reg(x,y,order='max')
#'
#' ## x-only paritioning (Univariate only)
#' NNS.reg(x,y,type="XONLY")
#'
#' ## Logistic Regression (Univariate only)
#' NNS.reg(x,y,type="LOGIT")
#'
#' ## For Multiple Regression:
#' x<-cbind(rnorm(100),rnorm(100),rnorm(100)); y<-rnorm(100)
#' NNS.reg(x,y,point.est=c(.25,.5,.75))
#'
#' ## For Multiple Regression based on Synthetic X* (Dimension Reduction):
#' x<-cbind(rnorm(100),rnorm(100),rnorm(100)); y<-rnorm(100)
#' NNS.reg(x,y,point.est=c(.25,.5,.75),type="CLASS")
#'
#' ## IRIS dataset example:
#' #Dimension Reduction:
#' NNS.reg(iris[,1:4],iris[,5],type="CLASS",order=5,dep.order=1)
#' #Multiple Regression:
#' NNS.reg(iris[,1:4],iris[,5],order=2)
#'
#' ## To call fitted values:
#' NNS.reg(x,y)$Fitted
#'
#' ## To call partial derivative (univariate regression only):
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.reg(x,y)$derivative
#'
#' @export


NNS.reg = function (x,y,
                    order=NULL,
                    s.t.n=.99,
                    type = NULL,
                    point.est = NULL,
                    location = 'top',
                    return.values = TRUE,
                    plot = TRUE,
                    residual.plot=TRUE,
                    threshold = 0,
                    dep.order=NULL,
                    precision="LOW",
                    n.best=2,
                    text=FALSE,
                    noise.reduction='mean',
                    norm=NULL,
                    dist="L2"){

  R2s = numeric()
  original.columns = ncol(x)
  original.variable = x
  np = nrow(point.est)

  if(!is.null(order)){
    if(order=='max'){
      order=ceiling(log2(length(y)))} else {order=order}} else {order=order}

  if(class(x)=='factor' | class(y)=='factor'){dep.order=1}

  if(!is.null(ncol(original.variable))){
    if(ncol(original.variable)==1){
      x=original.variable
    } else {
      if(is.null(type)){
        return(NNS.M.reg(x,y,point.est=point.est,plot=plot,residual.plot=plot,order=order,n.best=n.best,type=type,location=location,precision=precision,text=text,noise.reduction=noise.reduction,norm = norm,dist=dist))
      } # Multivariate NULL type
      else{
        if(type=="CLASS"){
          dep.order=dep.order
          x= data.matrix(x)
          y= as.numeric(y)

          x.star.dep = numeric()
          x.star.coef = numeric()
          x.star.matrix = matrix(nrow=length(y))

          for (i in 1:ncol(original.variable)){
            x.star.structure = NNS.dep(as.numeric(x[,i]),y,print.map = FALSE,order = dep.order)

            x.star.dep[i] = x.star.structure$Dependence
            x.star.coef[i]=  x.star.structure$Correlation

            if(abs(x.star.coef[i])<=threshold){
              x.star.coef[i]=0
            }

            x.star.matrix =  cbind(x.star.matrix,x.star.coef[i]*as.numeric(original.variable[,i]))

            if(i == ncol(original.variable)){
              synthetic.x.equation=(paste0("Synthetic Independent Variable X* = (",paste(format(x.star.coef[1:i],digits = 4),paste("X",1:i,sep = ''),sep='*',collapse = "  "),")/",sum(abs(x.star.coef)>0)))
            }
          } #ncol original variable

          if(!is.null(point.est)){
            new.point.est=numeric()
            if(is.null(np)){
              new.point.est= sum(point.est*x.star.coef)/sum(abs(x.star.coef)>0)}
            else {
              for(i in 1:np){
                new.point.est[i] = sum(point.est[i,]*x.star.coef)/sum(abs(x.star.coef)>0)
              }
              point.est=new.point.est
            }
          }

          x = rowSums(x.star.matrix[,2:(1+ncol(original.variable))])/ncol(original.variable)

        } # Multivariate "CLASS" type
      } #Multivariate Not NULL type
    }

  } #Multivariate

  if(is.null(original.columns)){
    synthetic.x.equation=NULL

    if(length(y)<100){dep.order=1}else{dep.order=dep.order}

    dependence = NNS.dep(x,y,print.map = FALSE)$Dependence

  } else {
    if(type=="CLASS") dependence=mean(x.star.dep)}

  if(is.null(order)){
    dep.reduced.order=ceiling(ceiling(log2(length(y)))*dependence)} else {
      dep.reduced.order=order
    }

    if(dependence>s.t.n ){
        if(is.null(type)){
            part.map = NNS.part(x,y,noise.reduction='off',order=dep.reduced.order)
        if(length(part.map$regression.points[,1])==0){
            part.map=NNS.part(x,y,noise.reduction='off',order = min(nchar(part.map$df$quadrant)),overfit = T)
        } else {part.map=part.map
          }
        } # NULL type

      if(!is.null(type)){
        part.map=NNS.part(x,y,type = "XONLY",noise.reduction='off',order = dep.reduced.order)

        if(length(part.map$regression.points[,1])==0){
          part.map=NNS.part(x,y,noise.reduction='off',type="XONLY",order = min(nchar(part.map$df$quadrant)),overfit = T)
        } else {part.map=part.map
        }
      } # type
    } # Dependence > s.t.n

    if(dependence<=s.t.n){
      if(is.null(type)){
        part.map = NNS.part(x,y,noise.reduction=noise.reduction,type = "XONLY", order=dep.reduced.order)

        if(length(part.map$regression.points[,1])==0){
          part.map=NNS.part(x,y,type = "XONLY",noise.reduction=noise.reduction,order = min(nchar(part.map$df$quadrant)),overfit = T)
        } else {part.map=part.map
        }
      } # NULL type


      if(!is.null(type)){
        part.map = NNS.part(x,y,type = "XONLY",noise.reduction=noise.reduction, order = dep.reduced.order)

        if(length(part.map$regression.points[,1])==0){
          part.map=NNS.part(x,y,type = "XONLY",noise.reduction=noise.reduction,order = min(nchar(part.map$df$quadrant)),overfit = T)
        } else {part.map=part.map
        }
      } # type
    } # Dependence < s.t.n

    naive.order =min(nchar(part.map$df$quadrant))-1

    Regression.Coefficients = data.frame(matrix(ncol=3))
    colnames(Regression.Coefficients) = c('Coefficient','X Lower Range','X Upper Range')

    regression.points=part.map$regression.points


    min.range = min(regression.points[,1])
    max.range = max(regression.points[,1])

    mode=function(x) {
      if(length(x)>1){
        d <- density(x)
        d$x[which.max(d$y)]
      }else{x}
    }

    Dynamic.average.min = mean(y[x<=min.range])
    Dynamic.average.max = mean(y[x>=max.range])

    ###Endpoints
    if(length(x[x<min.range])>0){
      if(dependence<s.t.n){
        x0 = Dynamic.average.min} else {
          x0 = unique(y[x==min(x)])} }  else {x0 = unique(y[x==min(x)])}

    if(length(x[x>max.range])>0){
      if(dependence<s.t.n){x.max = Dynamic.average.max} else {x.max = unique(y[x==max(x)])}}  else { x.max = unique(y[x==max(x)])}


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

    point.est.y=numeric()

    for (i in 1:p){

      z=(which(x>=Regression.Coefficients[i,2] & x<Regression.Coefficients[i,3]))

      z.diff = ((x[z]- Regression.Coefficients[i,2])*Regression.Coefficients[i,1])+regression.points[i,2]

      fitted.new =  cbind(z,z.diff)

      fitted = rbind(fitted,fitted.new)
      fitted = fitted[order(fitted[,1]),]

  ###Point estimates
      if(is.null(point.est)){point.est.y = NULL}
      else{
          for(j in 1:length(point.est)){
              if(point.est[j]>=Regression.Coefficients[i,2] & point.est[j]<Regression.Coefficients[i,3]){
                  point.est.y[j] = (point.est[j] - Regression.Coefficients[i,2])*(Regression.Coefficients[i,1])+regression.points[i,2]
              } else {
                    if(point.est[j]<Regression.Coefficients[1,2]){
                        point.est.y[j] = ((point.est[j] - Regression.Coefficients[1,2])*(Regression.Coefficients[1,1]))+(regression.points[1,2])
                    } else {
                          if(point.est[j]>Regression.Coefficients[p,2]){
                                point.est.y[j] = ((point.est[j] - Regression.Coefficients[(p-0),2])*(Regression.Coefficients[(p-1),1]))+(regression.points[(p),2])
                          }
                      }
                  }
          } #j loop
      } # else NULL point.est
    } #for i in p loop

    Values = (cbind(x,Fitted=fitted[,2],Actual=y,Difference=fitted[,2]-(y),
                    Accuracy=abs(round(fitted[,2])-(y))))

    MSE = mean((fitted[,2]-y)^2)
    y.fitted=fitted[,2]
    Fitted.values = y.fitted
    Prediction.Accuracy=(length(y)-sum(abs(round(y.fitted)-(y))>0))/length(y)

    R2=  (sum((fitted[,2]-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((fitted[,2]-mean(y))^2))

    R2.adj = R2

    ###Plotting and regression equation
    if(plot==TRUE){
      r2.leg=bquote(bold(R^2 == .(format(R2,digits=4))))
      xmin= min(c(point.est,x))
      xmax= max(c(point.est,x))
      ymin= min(c(point.est.y,y))
      ymax= max(c(point.est.y,y))

      plot(x,y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),col='steelblue',main=paste(paste0("NNS Order = ",naive.order)
                                                                              ,sep="\n"),
           xlab = if(!is.null(original.columns))
           {if(original.columns>1){paste0("Synthetic X* ","(Segments = ",(p-1),')')}}else{paste0("X  ","(Segments = ",(p-1),")",sep='')},
           ylab="Y",mgp=c(2.5,0.5,0),
           cex.lab=2,cex.main=2)

      ### Plot Regression points and fitted values and legend
      points(na.omit(regression.points[order(regression.points),]),col='red',pch=15)
      lines(na.omit(regression.points[order(regression.points),]),col='red',lwd=2,lty = 2)


      if(!is.null(point.est)){ points(point.est,point.est.y,
                                      col='green',pch=18,cex=1.5)

        legend(location, bty="n", y.intersp = 0.75,legend=r2.leg)}

      if(is.null(point.est)){
        legend(location, bty="n", y.intersp = 0.75,legend=r2.leg)
      }

      if(!is.null(point.est)){
        for(i in 1:length(point.est)){
          if(point.est[i]>max(x)) segments(point.est[i],point.est.y[i],regression.points[p+ceiling(abs(p-q)/2),1],regression.points[p+ceiling(abs(p-q)/2),2],col="green",lty=2)
          if(point.est[i]<min(x)) segments(point.est[i],point.est.y[i],regression.points[1,1],regression.points[1,2],col="green",lty=2)
        } }

    }# plot TRUE bracket

### Return Values
    if(return.values == TRUE){
        return(list("R2"=R2, "MSE"=MSE, "Prediction.Accuracy"=Prediction.Accuracy,"equation"=synthetic.x.equation,"Segments" = p-1, "derivative"=Regression.Coefficients[-p,],"Point"=point.est,"Point.est"=point.est.y,"regression.points"=regression.points,"partition"=cbind(part.map$df[,2:3],"Y.hat"=as.numeric(y.fitted)),"Fitted"=as.numeric(y.fitted),"Fitted.xy"=cbind(x,"Y.hat"=as.numeric(y.fitted))))
  }

}

#' NNS Multivariate Regression (INTERNAL CALL FOR \link{NNS.reg})
#'
#' Called by \code{NNS.reg} for multivariate regression analysis.
#'
#'
#' @param B Complete dataset of independent variables (IV) in matrix form.
#' @param y Dependent variable (DV).
#' @param order Controls the number of the \code{NNS.reg}.
#' @param s.t.n Signal to noise parameter, sets the threshold of \code{NNS.dep} which reduces \code{"order"} when \code{order=NULL}.  Defaults to 0.9 to ensure high dependence for higher \code{"order"} and endpoint determination.
#' @param n.best Sets the number of closest regression points to use in kernel weighting.  Defaults to 1.  Should be validated on hold-out set in conjunction with \code{"precision"} parameter.
#' @param type Controls the partitioning in \code{NNS.reg}.  Defaults to \code{type="XONLY"} for IV based partitioning.   \code{type=NULL} for both IV and DV partitioning.
#' @param point.est Generates a fitted value of \code{y} for a vector or matrix of IV coordinates.
#' @param plot Generates a 3d scatter plot with regression points using \link{plot3d}
#' @param residual.plot Generates a \code{matplot} for Y.hat and Y
#' @param location Sets the location of the legend
#' @param precision  Increases speed of computation at the expense of precision.  2 settings offered: \code{"LOW"} (Default setting), and \code{"HIGH"}.  \code{"HIGH"} is the limit condition of every observation as a regression point.
#' @param text If performing a text classification, set \code{text=TRUE}.  Defaults to FALSE.
#' @param noise.reduction In low signal:noise situations,\code{noise.reduction="mean"}  uses means for \link{NNS.dep} restricted partitions, \code{noise.reduction="median"} uses medians instead of means for \link{NNS.dep} restricted partitions, while \code{noise.reduction="mode"}  uses modes instead of means for \link{NNS.dep} restricted partitions.  \code{noise.reduction=NULL} (Default setting) allows for maximum possible fit and specific \code{order} specification.
#' @param norm Normalizes regressors between 0 and 1 for multivariate regression when set to \code{norm="std"}, or normalizes regressors according to \link{NNS.norm} when set to \code{norm="NNS"}. Defaults to NULL.
#' @return Returns the values: \code{"Fitted"} for only the fitted values of the DV; \code{"regression.points"} provides the points for each IV used in the regression equation for the given order of partitions; \code{"rhs.partitions"} returns the partition points for each IV; \code{"partition"} returns the DV, quadrant assigned to the observation and fitted value, and  \code{"Point.est"} for predicted values.
#' @keywords  multiple nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}


NNS.M.reg <- function (B,y,order=NULL,s.t.n=0.9,n.best=1,type=NULL,point.est=NULL, plot=FALSE,residual.plot=TRUE,location=NULL,precision="LOW",text=FALSE,noise.reduction=FALSE,norm=NULL){

  if(is.null(ncol(B))){B=t(t(B))}
  n=ncol(B)

  np=nrow(point.est)
  if(is.null(np)&!is.null(point.est)){point.est=t(point.est)
  }else{point.est=point.est}
  findInterval2Vec = Vectorize(findInterval2,vectorize.args="x")
  ### For Multiple regressions
  ###  Turn each column into numeric values
  original.variable = B

  new.variable=data.matrix(original.variable)
  B=new.variable
  if(!is.null(norm)){
  if(norm=='std'){
  B=apply(B,2,function(b) (b-min(b))/(max(b)-min(b)))}else{
  B=NNS.norm(B)}
  if(!is.null(point.est)){
  if(norm=='std'){
  point.est=apply(point.est,2,function(c) (c-min(rbind(c,original.variable)))/(max(rbind(c,original.variable))-min(rbind(c,original.variable))))}else{
  point.est=NNS.norm(rbind(point.est,original.variable))[1:np,]}
  }}else{B=B
  point.est=point.est
  }

  colnames(B)=colnames(original.variable)
  y=as.numeric(y)


  reg.points=list()
  sections = list()


  ###  Regression Point Matrix

  for(i in 1:n){
      reg.points[[i]] = NNS.reg(B[,i],y,order=order,type=type,noise.reduction=NULL,precision = precision,plot = F)$regression.points[,1]
  }

  reg.points.matrix=do.call('cbind',lapply(reg.points, `length<-`, max(lengths(reg.points))))

  if(length(reg.points.matrix[,1])==0){
    for(i in 1:n){
      part.map=partition.map(B[,i],y,order=order,type=type,noise.reduction=noise.reduction)
      dep=NNS.dep(B[,i],y,order=1)$Dependence
    if(dep>s.t.n){
        reg.points[[i]] = partition.map(B[,i],y,order=round(dep*max(nchar(part.map$df[,3]))),type=type,noise.reduction=NULL)$regression.points[,1]}
      else{reg.points[[i]] = partition.map(B[,i],y,order=round(dep*max(nchar(part.map$df[,3]))),noise.reduction=noise.reduction,type="XONLY")$regression.points[,1]}
    }


  reg.points.matrix=do.call('cbind',lapply(reg.points, `length<-`, max(lengths(reg.points))))}


  if(is.null(colnames(B))){
    colnames.list=list()
    for(i in 1:n){
      colnames.list[i]=paste0("X",i)
    }
    colnames(reg.points.matrix)=c(colnames.list)}
  else{
    colnames(reg.points.matrix)=colnames(B)}
  reg.length=length(na.omit(reg.points.matrix[,1]))

  ### Find intervals in regression points for each variable
  intervals = list()

  for(j in 1:n){
    intervals[[j]]=findInterval2Vec(B[,j],sort(unique(reg.points.matrix[,j])))
  }

  intervals = matrix(unlist(intervals),nrow = length(y),ncol = n)


  ## Create unique identifier of each observation's interval
  intervals = apply(intervals, 1 , paste , collapse = "" )
  y.identifier = cbind(as.numeric(y),as.numeric(intervals))
  if(text==TRUE){
   if(n<10){y.identifier = cbind(as.numeric(y),as.numeric(intervals))}else
     {y.identifier = cbind(as.numeric(y),as.numeric(y))}
  }

  ## Match y to unique identifier
  overall.matrix = cbind(B,y,y.identifier[,2])
  overall.matrix = overall.matrix[order(overall.matrix[,(n+2)]),]
  y.hat.matrix =  data.frame(matrix(ncol = 2))
  y.hat=numeric()

  for(i in 1:length(y)){
    y.hat[i] = mean(c(y.identifier[,1][y.identifier[,2]==y.identifier[i,2]]))
  }

  y.identifier = cbind(y.identifier,y.hat)
  y.identifier = y.identifier[order(y.identifier[,2]),]

  colnames(y.identifier)=c('y','NNS identifier','fitted value')

  fitted.matrix = cbind(B,y.hat)

  unique.identifiers=unique(overall.matrix[,(n+2)])

  ### GET REGRESSION POINTS USING Y.IDENTIFIER!!!
  REGRESSION.POINT.MATRIX = data.frame(matrix(ncol=n+2))
  colnames(REGRESSION.POINT.MATRIX)=c(colnames(reg.points.matrix),"Y")
  for(j in 1:length(unique.identifiers)){
    for(i in 1:n){
      REGRESSION.POINT.MATRIX[j,i]=mean(c(overall.matrix[overall.matrix[,(n+2)]==unique.identifiers[j],i]))
      REGRESSION.POINT.MATRIX[j,(n+1)]=mean(c(overall.matrix[overall.matrix[,(n+2)]==unique.identifiers[j],(n+1)]))
      REGRESSION.POINT.MATRIX[j,(n+2)]=mean(c(overall.matrix[overall.matrix[,(n+2)]==unique.identifiers[j],(n+2)]))
    }
  }
  REGRESSION.POINT.MATRIX=na.omit(REGRESSION.POINT.MATRIX)


  ### DISTANCES
  ## Calculate distance from each point in REGRESSION.POINT.MATRIX
  if(!is.null(point.est)){
    distance<- function(dist.est){
      distances=sweep(REGRESSION.POINT.MATRIX[,(1:n)],2,as.numeric(dist.est))
      distances[distances==0]<- 1e-10
      row.sums=as.numeric(sqrt(rowSums(distances^2)))
      distances.max=sweep(overall.matrix[,(1:n)],2,as.numeric(dist.est))
      distances.max[distances.max==0]<- 1e-10
      row.sums.max=as.numeric(sqrt(rowSums(distances.max^2)))

      total.row.sums = sum(1/row.sums)
      total.row.sums.max = sum(1/row.sums.max)

      weights = (1/row.sums)/total.row.sums
      weights.max.prec = (1/row.sums.max)/total.row.sums.max

      which_nth_highest <- function(xx, n)
      {
        ux <- unique(xx)
        nux <- length(ux)
        which(xx == sort(ux, partial = nux - n + 1)[nux - n + 1])
      }

      if(precision=="LOW"){
          highest=numeric()

          for(k in 1:n.best){
              highest[k]=which_nth_highest(weights,k)[1]
          }

          weights[-c(highest)]<-0
          weights.sum=sum(weights)

          for(l in 1:length(weights)){
              weights[l]=weights[l]/weights.sum
          }

          single.estimate = sum(weights*REGRESSION.POINT.MATRIX[,(n+1)])
      }


      if(precision=="HIGH"){
          highest=numeric()
          for(k in 1:n.best){
          highest[k]=which_nth_highest(weights.max.prec,k)[1]
        }
        weights.max.prec[-c(highest)]<-0
        weights.max.sum=sum(weights.max.prec)


        for(l in 1:length(weights.max.prec)){
          weights.max.prec[l]=weights.max.prec[l]/weights.max.sum
        }

        single.estimate = sum(weights.max.prec*overall.matrix[,(n+1)])

      }


     return(single.estimate)
    }


  ### Point estimates

    predict.fit=numeric()
    predict.fit.iter=numeric()
    if(is.null(np)){
      predict.fit = distance(dist.est = point.est)
    }
    if(!is.null(np)){
      for(j in 1:np){
        predict.fit.iter[j] = distance(dist.est = as.vector(point.est[j,]))
      }
      predict.fit=predict.fit.iter
    }

 }#is.null point.est


  R2=  (sum((y.hat-mean(y))*(y-mean(y)))^2)/(sum((y-mean(y))^2)*sum((y.hat-mean(y))^2))



  if(plot==TRUE&&n==2){


    REGRESSION.POINT.MATRIX=REGRESSION.POINT.MATRIX[order(REGRESSION.POINT.MATRIX[,(n+1)]),]

    plot3d(x=B[,1],y=B[,2],z=y,box=F,size = 3,col='steelblue',xlab="X1", ylab="X2", zlab="Y" )
    points3d(x=REGRESSION.POINT.MATRIX[,1],y=REGRESSION.POINT.MATRIX[,2],z=REGRESSION.POINT.MATRIX[,3],col='red',size=5)


    if(!is.null(point.est)){
      if(is.null(np)){
        points3d(x=point.est[1],y=point.est[2],z=predict.fit,col='green',size=5)
      } else {
        points3d(x=point.est[,1],y=point.est[,2],z=predict.fit,col='green',size=5)
      }
    }


  }

  if(residual.plot==TRUE){
    r2.leg=bquote(bold(R^2 == .(format(R2,digits=4))))
    matplot(cbind(y,y.hat),type = 'l',xlab="Index",ylab="Y (black) and Y.hat (red)",cex.lab=1.5,mgp=c(2,.5,0))

    title(main = paste0("NNS Order = ",order),cex.main=2)
    legend(location,legend =r2.leg,bty = 'n')

  }

  if(n==2){
    if(!is.null(point.est)){
      return(list(R2=R2,Fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier,Point.est=predict.fit,Fitted.xy=fitted.matrix))} else {return(list(R2=R2,Fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier,Fitted.xy=fitted.matrix))}
  }

  else{if(!is.null(point.est)){
    return(list(R2=R2,Fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier,Point.est=predict.fit,Fitted.xy=fitted.matrix))} else {return(list(R2=R2,Fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))],partition=y.identifier,Fitted.xy=fitted.matrix))}
  }



}

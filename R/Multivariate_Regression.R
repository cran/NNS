#' VN Multivariate Regression (INTERNAL CALL FOR \link{VN.reg})
#'
#' Called by \code{VN.reg} for multivariate regression analysis.
#'
#'
#' @param B Complete dataset of independent variables (IV) in matrix form.
#' @param y Dependent variable (DV).
#' @param order Controls the number of the \code{VN.reg}.  Defaults to \code{order=1}.
#' @param s.t.n Signal to noise parameter, sets the threshold of \code{VN.dep} which reduces \code{"order"} when \code{order=NULL}.  Defaults to 0.9 to ensure high dependence for higher \code{"order"} and endpoint determination.
#' @param n.best Sets the number of closest regression points to use in kernel weighting.  Defaults to number of independent variables.
#' @param type Controls the partitioning in \code{VN.reg}.  Defaults to \code{type="XONLY"} for IV based partitioning.   \code{type=NULL} for both IV and DV partitioning.
#' @param point.est Generates a fitted value of \code{y} for a vector or matrix of IV coordinates.
#' @param plot Generates a 3d scatter plot with regression points using \link{plot3d}
#' @param residual.plot Generates a \code{matplot} for Y.hat and Y
#' @param location Sets the location of the legend
#' @return Returns the values: \code{"Fitted"} for only the fitted values of the DV; \code{"regression.points"} provides the points for each IV used in the regression equation for the given order of partitions; \code{"rhs.partitions"} returns the partition points for each IV; \code{"partition"} returns the DV, quadrant assigned to the observation and fitted value, and  \code{"Point.est"} for predicted values.
#' @keywords  multiple nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}


VN.M.reg <- function (B,y,order=1,s.t.n=0.9,n.best=NULL,type="XONLY",point.est=NULL, plot=FALSE,residual.plot=TRUE,location=location){

  if(is.null(ncol(B))){B=t(t(B))}
  n=ncol(B)

  if(is.null(n.best)){n.best=n}else{n.best=n.best}
  np=nrow(point.est)

  findInterval2Vec = Vectorize(findInterval2,vectorize.args="x")
  ### For Multiple regressions

  original.variable = B
  new.variable = matrix(nrow=nrow(B))

  ###  Turn each column into numeric values
  for (i in 1:ncol(original.variable)){
    new.variable = cbind(new.variable,as.numeric(original.variable[,i]))
  }

  B=new.variable[,-1]
  colnames(B)=colnames(original.variable)
  y=as.numeric(y)


  reg.points=list()
  sections = list()


  ###  Regression Point Matrix
  reg.points.matrix = matrix(ncol=n)
  for(i in 1:n){
    reg.points[[i]] = VN.reg(B[,i],y,plot = FALSE,order=order,type = "XONLY",s.t.n=s.t.n)$regression.points[,1]
  }

  reg.points.matrix=do.call('cbind',reg.points)
  reg.length=length(reg.points.matrix[,1])

  ### Find intervals in regression points for each variable
  intervals = list()

  for(j in 1:n){
    intervals[[j]]=findInterval2Vec(B[,j],sort(unique(reg.points.matrix[,j])))
  }

  intervals = matrix(unlist(intervals),nrow = length(y),ncol = n)


  ## Create unique identifier of each observation's interval
  intervals = apply(intervals, 1 , paste , collapse = "" )
  y.identifier = cbind(as.numeric(y),as.numeric(intervals))     ## Match y to unique identifier
  overall.matrix = cbind(B,y,y.identifier[,2])

  y.hat.matrix =  data.frame(matrix(ncol = 2))
  y.hat=numeric()

  for(i in 1:length(y)){
    y.hat[i] = mean(c(y.identifier[,1][y.identifier[,2]==y.identifier[i,2]]))
  }



  y.identifier = cbind(y.identifier,y.hat)

  y.identifier = y.identifier[order(y.identifier[,2]),]

  colnames(y.identifier)=c('y','quadrant','fitted value')

  fitted.matrix = cbind(B,y.hat)

  unique.identifiers=unique(overall.matrix[,(n+2)])

  ### GET REGRESSION POINTS USING Y.IDENTIFIER!!!
  REGRESSION.POINT.MATRIX = data.frame(matrix(ncol=n+2))
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
      distances=data.frame(nrow=length(REGRESSION.POINT.MATRIX[,1]),ncol=n)
      row.sums=numeric()

      for(i in 1:length(REGRESSION.POINT.MATRIX[,1])){
        for(j in 1:n){
          distances[i,j]=abs(dist.est[j]-REGRESSION.POINT.MATRIX[i,j])
        }

        distances[i,][distances[i,]==0]<- 1e-10
        row.sums[i] = sum(distances[i,]^2)
      }

      total.row.sums = sum(1/row.sums)

      weights = (1/row.sums)/total.row.sums

      which_nth_highest <- function(xx, n)
      {
        ux <- unique(xx)
        nux <- length(ux)
        which(xx == sort(ux, partial = nux - n + 1)[nux - n + 1])
      }


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
      return(single.estimate)
    }
  }

  ### Point estimates
  if(!is.null(point.est)){
    predict.fit=numeric()
    predict.fit.iter=numeric()
    if(is.null(np)){
      predict.fit = distance(dist.est = point.est)
    }
    if(!is.null(np)){
      for(j in 1:np){
        predict.fit.iter[j] = distance(dist.est = point.est[j,])
      }
      predict.fit=predict.fit.iter
    }

  }


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

    # movie3d(spin3d(axis = c(.75, .75, .75)), duration = 20,dir = getwd(),convert=FALSE)

  }

  if(residual.plot==TRUE){
    r2.leg=bquote(bold(R^2 == .(format(R2,digits=4))))
    matplot(cbind(y,y.hat),type = 'l',xlab="Index",ylab="Y (black) and Y.hat (red)",cex.lab=1.5,mgp=c(2,.5,0))

    title(main = paste0("NNS Order = ",order),cex.main=2)
    legend('top',legend =r2.leg,bty = 'n')

  }

  if(n==2){
    if(!is.null(point.est)){
      return(list(R2=R2,Fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier,Point.est=predict.fit))} else {return(list(R2=R2,Fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier))}
  }

  else{if(!is.null(point.est)){
    return(list(R2=R2,Fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier,Point.est=predict.fit))} else {return(list(R2=R2,Fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier))}
  }



}

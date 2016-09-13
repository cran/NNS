#' VN Multivariate Regression
#'
#' Called by \code{VN.reg} for multivariate regression analysis.
#'
#'
#' @param B Complete dataset of independent variables (IV) in matrix form.
#' @param y Dependent variable (DV).
#' @param order Controls the number of the \code{VN.reg}.  Defaults to \code{order=1}.
#' @param type Controls the partitioning in \code{VN.reg}.  Defaults to \code{type="XONLY"} for IV based partitioning.   \code{type=NULL} for both IV and DV partitioning.
#' @param point.est Generates a fitted value of \code{y} for a vector or matrix of IV coordinates.
#' @param plot Generates a 3d scatter plot with regression points using \link{plot3d}
#' @param residual.plot Generates a \code{matplot} for Y.hat and Y
#' @keywords  multiple nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123);size=30;x_1=rnorm(size,mean=5,sd=1);x_2=rnorm(size,mean=1,sd=5); x_3=runif(size)
#' y=x_1^2+x_2+x_3
#' B=cbind(x_1,x_2,x_3)
#' VN.M.reg(B,y)
#' @export

VN.M.reg <- function (B,y,order=1,type="XONLY",point.est=NULL, plot=FALSE,residual.plot=TRUE){

  n=ncol(B)
  np=nrow(point.est)
  overall.matrix = cbind(B,y)

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
    reg.points[[i]] = VN.reg(B[,i],y,return.equation = TRUE,plot = FALSE,order=order,type = "XONLY")$regression.points[,1]
  }

  reg.points.matrix=do.call('cbind',reg.points)
  reg.length=length(reg.points.matrix[,1])

  ### Find intervals in regression points for each variable
  intervals = list()

  for(j in 1:n){
    intervals[[j]]=findInterval2Vec(B[,j],reg.points.matrix[,j])
  }

  intervals = matrix(unlist(intervals),nrow = length(y),ncol = n)


  ## Create unique identifier of each observation's interval
  intervals = apply(intervals, 1 , paste , collapse = "" )
  y.identifier = cbind(as.numeric(y),as.numeric(intervals))     ## Match y to unique identifier
  y.hat.matrix =  data.frame(matrix(ncol = 2))
  y.hat=numeric()

  for(i in 1:length(y)){
    y.hat[i] = mean(c(y.identifier[,1][y.identifier[,2]==y.identifier[i,2]]))
  }


  y.identifier = cbind(y.identifier,y.hat)

  y.identifier = y.identifier[order(y.identifier[,2]),]

  colnames(y.identifier)=c('y','quadrant','fitted value')

  fitted.matrix = cbind(B,y.hat)

  ### Matrix is means of intervals, thus regression points.
  mult.reg.points.matrix = as.data.frame(matrix(ncol=n))

  # X1 and X2 are for plot3d
  combo.x1 = numeric()
  combo.x2 = numeric()
  for(j in 1:n){
    for(i in 1:(reg.length-1)){
      mult.reg.points.matrix[1,j]= mean(reg.points.matrix[1,j])
      mult.reg.points.matrix[(i+1),j]= mean(reg.points.matrix[(i:(i+1)),j])
      combo.x1[i]= mean(reg.points.matrix[(i:(i+1)),1])
      combo.x2[i]= mean(reg.points.matrix[(i:(i+1)),2])
      mult.reg.points.matrix[reg.length,j]= mean(reg.points.matrix[reg.length,j])
    }
  }

  ### Expand the points to see all intervals
  e.g.mult.reg.points.matrix=expand.grid(mult.reg.points.matrix)

  ### Find y for unique intervals in Multiple Regression Points Matrix m.r.p.
  y.m.r.p.identifier=list()
  for(j in 1:n){
    y.m.r.p.identifier[[j]]=findInterval2Vec(e.g.mult.reg.points.matrix[,j],reg.points.matrix[,j])
  }

  y.m.r.p.identifier = matrix(unlist(y.m.r.p.identifier),nrow =length(e.g.mult.reg.points.matrix[,1]),ncol = n)
  y.m.r.p.identifier = apply(y.m.r.p.identifier, 1 , paste , collapse = "" )
  y.m.r.p.identifier = as.numeric(y.m.r.p.identifier)

  e.g.mult.reg.points.matrix = cbind(e.g.mult.reg.points.matrix, y.m.r.p.identifier)


  y.m.r.p.identifier = e.g.mult.reg.points.matrix[e.g.mult.reg.points.matrix[,(n+1)]%in%y.identifier[,2],(n+1)]

  e.g.mult.reg.points.matrix = e.g.mult.reg.points.matrix[e.g.mult.reg.points.matrix[,(n+1)]%in%y.identifier[,2],-(n+1)]

  Conditional.y.means = numeric()

  for(i in 1:length(y.m.r.p.identifier)){
    Conditional.y.means[i] = mean(y.identifier[,1][y.identifier[,2]==y.m.r.p.identifier[i]])
  }

  REGRESSION.POINT.MATRIX=cbind(e.g.mult.reg.points.matrix, Conditional.y.means, y.m.r.p.identifier)


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
        row.sums[i] = prod(distances[i,]^2)

      }


      total.row.sums = sum(1/row.sums)

      weights = (1/row.sums)/total.row.sums

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
    # For 3d plot which can only have 2 RHS variables
    max.x1 = max(overall.matrix[,1])
    max.x2 = max(overall.matrix[,2])
    REGRESSION.POINT.MATRIX=rbind(REGRESSION.POINT.MATRIX,c(max.x1,overall.matrix[,2][overall.matrix[,1]==max.x1],overall.matrix[,(n+1)][overall.matrix[,1]==max.x1]))
    REGRESSION.POINT.MATRIX=rbind(REGRESSION.POINT.MATRIX,c(overall.matrix[,1][overall.matrix[,2]==max.x2],max.x2,overall.matrix[,(n+1)][overall.matrix[,2]==max.x2]))
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
    matplot(cbind(y,y.hat),type = 'l',ylab="Y and Y.hat")
    legend('top',legend = c("Y (black)","NNS Fit (red)",paste("R2= ",format(R2,digits = 6),sep = '')),col=c('black','red'),bty = 'n')
  }

  if(n==2){
    if(!is.null(point.est)){
      return(list(R2=R2,fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier,prediction=predict.fit))} else {return(list(R2=R2,fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier))}
  }

  else{if(!is.null(point.est)){
    return(list(R2=R2,fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier,prediction=predict.fit))} else {return(list(R2=R2,fitted=y.hat,rhs.partitions=reg.points.matrix, regression.points=REGRESSION.POINT.MATRIX[,(1:(n+1))] ,partition=y.identifier))}
  }



}

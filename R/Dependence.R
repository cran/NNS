#' NNS Dependence
#'
#' Returns the dependence between two variables based on higher order partial moment correlations measured by frequency or area.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param order Controls the level of quadrant partitioning.  Default to \code{order=NULL}.  Errors can generally be rectified by setting \code{order=1}.
#' @param degree Defaults to NULL to allow number of observations to be \code{degree} determinant.
#' @param print.map  Plots quadrant means.  Defaults to FALSE.
#' @return Returns the bi-variate \code{"Correlation"} and \code{"Dependence"} or correlation / dependence matrix for matrix input.
#' @keywords dependence, correlation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.dep(x,y)
#'
#' ## Correlation / Dependence Matrix
#' x<-rnorm(100); y<-rnorm(100); z<-rnorm(100)
#' B<-cbind(x,y,z)
#' NNS.cor(B)
#' @export

NNS.dep = function( x, y,order = NULL,
                   degree=NULL,
                   print.map=FALSE){


  if(is.null(degree)){degree=ifelse(length(x)<100,0,1)}else{degree=degree}
  if(!missing(y)){
      if(print.map==T){
          part.map = NNS.part(x,y,order=order, Voronoi=T)}
          else {
          part.map = NNS.part(x,y,order=order)
          }

  partitioned_df = part.map$df
  reg.points = part.map$regression.points

  clpm = numeric()
  cupm = numeric()
  dlpm = numeric()
  dupm = numeric()
  cor.rhos = numeric()
  dep.rhos = numeric()
  nonlin.cor = numeric()


  max.part = min(nchar(partitioned_df$quadrant))
  max.actual = max(nchar(partitioned_df$quadrant))
  part = nchar(partitioned_df$quadrant)

  if(max.part==max.actual){reduction=max.part-1}else{reduction=max.part}

  prior.partitioned_df=partitioned_df

  prior.partitioned_df[,'quadrant']=substr(partitioned_df$quadrant,1,reduction)

  int.points=numeric()

  for(item in unique(prior.partitioned_df$quadrant)){

    sub_prior=prior.partitioned_df[prior.partitioned_df$quadrant == item,]
    sub_x=sub_prior$x
    sub_y=sub_prior$y
    if(length(sub_x)==1){
        int.points=c(int.points,1)/length(x)
        dep.weight=0
        } else {
            dep.weight=length(sub_x)/length(x)
          }

      sub.clpm=Co.LPM(degree,degree,sub_x,sub_y)
      sub.cupm=Co.UPM(degree,degree,sub_x,sub_y)
      sub.dlpm=D.LPM(degree,degree,sub_x,sub_y)
      sub.dupm=D.UPM(degree,degree,sub_x,sub_y)

      clpm = c(clpm, sub.clpm)
      cupm = c(cupm, sub.cupm)
      dlpm = c(dlpm, sub.dlpm)
      dupm = c(dupm, sub.dupm)

      dep.rhos =  c(dep.rhos,dep.weight*abs((sub.clpm+sub.cupm-sub.dlpm-sub.dupm) / (sub.clpm+sub.cupm+sub.dlpm+sub.dupm)))

      nonlin.cor = c(nonlin.cor,dep.weight*(sub.clpm+sub.cupm-sub.dlpm-sub.dupm) / (sub.clpm+sub.cupm+sub.dlpm+sub.dupm))

}


  nonlin.cor = sum(na.omit(nonlin.cor))
  if(nonlin.cor<0) nonlin.cor=nonlin.cor-sum(int.points) else nonlin.cor=nonlin.cor+sum(int.points)

  if(is.na(nonlin.cor)){nonlin.cor=0}

  dep = sum(na.omit(dep.rhos))+sum(int.points)

  if(is.na(dep)){dep=0}


  return(list("Correlation"=nonlin.cor,"Dependence"= dep ))

  }#Not missing Y

  else{
  NNS.dep.matrix(x)
  }

}



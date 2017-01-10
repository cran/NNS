#' NNS Dependence
#'
#' Returns the dependence between two variables based on higher order partial moment correlations measured by frequency or area.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param order Controls the level of quadrant partitioning.  Default to \code{order=2}.  Errors can generally be rectified by setting \code{order=1}.
#' @param degree Defaults to NULL to allow number of observations to be \code{degree} determinant.
#' @param print.map  Plots quadrant means.  Defaults to FALSE.
#' @return Returns the \code{"Correlation"} and \code{"Dependence"}
#' @keywords dependence, correlation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.dep(x,y)
#' @export

NNS.dep = function( x, y,order = 2,
                   degree=NULL,
                   print.map=FALSE){


  if(is.null(degree)){degree=ifelse(length(x)<100,0,1)}else{degree=degree}

  if(print.map==T){
  part.map = partition.map(x,y,order=order, Voronoi=T)} else {
    part.map = partition.map(x,y,order=order)
  }

  partitioned_df = part.map$df
  reg.points = part.map$regression.points

  clpm = numeric(0)
  cupm = numeric(0)
  dlpm = numeric(0)
  dupm = numeric(0)
  cor.rhos = numeric(0)
  dep.rhos = numeric(0)
  nonlin.cor = numeric(0)


  max.part = min(nchar(partitioned_df$quadrant))
  max.actual = max(nchar(partitioned_df$quadrant))
  part = nchar(partitioned_df$quadrant)

  if(max.part==max.actual){reduction=max.part-1}else{reduction=max.part}

  prior.partitioned_df=partitioned_df

  prior.partitioned_df[,'quadrant']=substr(partitioned_df$quadrant,1,reduction)


  for(item in unique(prior.partitioned_df$quadrant)){

    sub_x=prior.partitioned_df[prior.partitioned_df$quadrant == item,'x']
    sub_y=prior.partitioned_df[prior.partitioned_df$quadrant == item,'y']

    weight=length(sub_x)/length(x)

      clpm = c(clpm, Co.LPM(degree,sub_x,sub_y))
      cupm = c(cupm, Co.UPM(degree,sub_x,sub_y))
      dlpm = c(dlpm, D.LPM(degree,degree,sub_x,sub_y))
      dupm = c(dupm, D.UPM(degree,degree,sub_x,sub_y))

      dep.rhos =  c(dep.rhos,weight*abs((Co.LPM(degree,sub_x,sub_y)+Co.UPM(degree,sub_x,sub_y)-D.LPM(degree,degree,sub_x,sub_y)-D.UPM(degree,degree,sub_x,sub_y)) / (Co.LPM(degree,sub_x,sub_y)+Co.UPM(degree,sub_x,sub_y)+D.LPM(degree,degree,sub_x,sub_y)+D.UPM(degree,degree,sub_x,sub_y))))

      nonlin.cor = c(nonlin.cor,weight*(Co.LPM(degree,sub_x,sub_y)+Co.UPM(degree,sub_x,sub_y)-D.LPM(degree,degree,sub_x,sub_y)-D.UPM(degree,degree,sub_x,sub_y)) / (Co.LPM(degree,sub_x,sub_y)+Co.UPM(degree,sub_x,sub_y)+D.LPM(degree,degree,sub_x,sub_y)+D.UPM(degree,degree,sub_x,sub_y)))

}

  nonlin.cor = sum(na.omit(nonlin.cor))

  if(is.na(nonlin.cor)){nonlin.cor=0}

  dep = sum(na.omit(dep.rhos))

  if(is.na(dep)){dep=0}


  return(list("Correlation"=nonlin.cor,"Dependence"= dep ))

}

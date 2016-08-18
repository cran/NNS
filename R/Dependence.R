#' VN Dependence
#'
#' Returns the dependence between two variables based on higher order partial moment correlations measured by frequency or area.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param order Controls the level of quadrant partitioning.  Defualts to NULL to have data determine maximum partitions.  Setting a specific order will activate \code{overide=TRUE} from \link{partition.map}
#' @param degree Defaults to NULL to allow number of observations to be \code{degree} determinant.
#' @param print.map  Displays partition mapping onto plot.  Defaults to TRUE.
#' @keywords dependence, correlation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' VN.dep(x,y)
#' @export

VN.dep = function( x, y,order = NULL,
                   degree=NULL,
                   print.map=TRUE){

  if(is.null(degree)){degree=ifelse(length(x)<100,0,1)}else{degree=degree}

  part.map = partition.map(x,y,order=order)

  partitioned_df = part.map$df
  reg.points = part.map$regression.points

  clpm = numeric(0)
  cupm = numeric(0)
  dlpm = numeric(0)
  dupm = numeric(0)
  cor.rhos = numeric(0)
  dep.rhos = numeric(0)

  if(print.map==TRUE){
    plot(x,y,col='blue',pch=20)
    abline(h=mean(y),v=mean(x),lwd=3,col='azure4')
    points(reg.points[,1],reg.points[,2],pch=19,lwd=2,col='red')}


  max.part = min(nchar(partitioned_df$master_part))
  part = nchar(partitioned_df$master_part)

  partitioned_df[,'master_part']=ifelse(part>max.part,
    substr(partitioned_df$master_part[part>max.part],1,nchar(partitioned_df$master_part[part>max.part])-1),partitioned_df$master_part)

  prior.partitioned_df=partitioned_df

  prior.partitioned_df[,'master_part']=substr(prior.partitioned_df$master_part,1,nchar(prior.partitioned_df$master_part)-1)


  for(item in unique(prior.partitioned_df$master_part)){

      sub_x = prior.partitioned_df[prior.partitioned_df$master_part == item, 'x']
      sub_y = prior.partitioned_df[prior.partitioned_df$master_part == item, 'y']


      clpm = c(clpm, Co.LPM(degree, mean(sub_x),mean(sub_y),sub_x,sub_y))
      cupm = c(cupm, Co.UPM(degree, mean(sub_x),mean(sub_y),sub_x,sub_y))
      dlpm = c(dlpm, D.LPM(degree,degree, mean(sub_x),mean(sub_y),sub_x,sub_y))
      dupm = c(dupm, D.UPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))


  }##nchar

  dep.rhos =  abs((clpm+cupm-dlpm-dupm) / (clpm+cupm+dlpm+dupm))

  m<- rbind((sum(clpm) +sum(cupm) -sum(dlpm) -sum(dupm))/(sum(clpm)+sum(cupm)+sum(dlpm)+sum(dupm)),

            sum(na.omit(dep.rhos))/length(na.omit(dep.rhos)))

  rownames(m) = c("Correlation","Dependence")


  m

}

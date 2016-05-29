#' VN Dependence
#'
#' Returns the dependence between two variables based on higher order partial moment correlations measured by frequency or area.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param order Controls the level of quadrant partitioning.  Defualts to NULL, but lower levels should be called for large (n).
#' @param degree Defaults to 0 for smaller number of observations.
#' @param print.map  Displays partition mapping onto plot.  Defaults to TRUE.
#' @keywords dependence
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.dep(x,y)}
#' @export

VN.dep = function( x, y,order = NULL,
                   degree=0,
                   print.map=TRUE){

  if(is.null(order)){
    for (i in 2:floor(log(length(x),4))){

      if(min(nchar(partition.map(x,y,i)$master_part))==i)
        order=i-1
    }}


  partitioned_df = partition.map(x, y,order)



  clpm = numeric(0)
  cupm = numeric(0)
  dlpm = numeric(0)
  dupm = numeric(0)
  cor.rhos = numeric(0)
  dep.rhos = numeric(0)

  if(print.map==TRUE){
    plot(x,y,col='blue',pch=20)
    abline(h=mean(y),v=mean(x),lwd=3,col='azure4')}

  prior.partitioned_df = partitioned_df
  prior.partitioned_df[,'master_part'] = substr(partitioned_df$master_part, 1, nchar(partitioned_df$master_part)-1)

  partition.lengths = numeric()


  for(item in unique(prior.partitioned_df$master_part)){
    partition.lengths[item] = nchar(item)
  }



  for(item in unique(prior.partitioned_df$master_part)){

    if(nchar(item) == min(partition.lengths) && min(partition.lengths)==order+0){

      sub_x = prior.partitioned_df[prior.partitioned_df$master_part == item, 'x']
      sub_y = prior.partitioned_df[prior.partitioned_df$master_part == item, 'y']


      clpm = c(clpm, Co.LPM(degree, mean(sub_x),mean(sub_y),sub_x,sub_y))
      cupm = c(cupm, Co.UPM(degree, mean(sub_x),mean(sub_y),sub_x,sub_y))
      dlpm = c(dlpm, D.LPM(degree,degree, mean(sub_x),mean(sub_y),sub_x,sub_y))
      dupm = c(dupm, D.UPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))



      if(print.map==TRUE){

        if(mean(sub_x)<mean(x) && mean(sub_y)<mean(y) ){
          segments(mean(sub_x),max(sub_y),mean(sub_x),min(sub_y),lty=3,lwd=2,col='red')
          segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3,lwd=2,col='red')
        }

        if(mean(sub_x)<mean(x) && mean(sub_y)>mean(y)){
          segments(mean(sub_x),max(sub_y),mean(sub_x),min(sub_y),lty=3,lwd=2,col='red')
          segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3,lwd=2,col='red')
        }
        if(mean(sub_x)>mean(x) && mean(sub_y)<mean(y)){
          segments(mean(sub_x),max(sub_y),mean(sub_x),min(sub_y),lty=3,lwd=2,col='red')
          segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3,lwd=2,col='red')
        }
        if(mean(sub_x)>mean(x) && mean(sub_y)>mean(y)){
          segments(mean(sub_x),max(sub_y),mean(sub_x),min(sub_y),lty=3,lwd=2,col='red')
          segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3,lwd=2,col='red')
        }}

    }

  }##nchar



  for(i in 1:(4^order)){

    dep.rhos[i] =  abs((clpm[i]+cupm[i]-dlpm[i]-dupm[i]) / (clpm[i]+cupm[i]+dlpm[i]+dupm[i]))
  }

  m<- rbind((sum(clpm) +sum(cupm) -sum(dlpm) -sum(dupm))/(sum(clpm)+sum(cupm)+sum(dlpm)+sum(dupm)),

            sum(na.omit(dep.rhos))/length(na.omit(dep.rhos)))

  rownames(m) = c("Correlation","Dependence")


  m



}

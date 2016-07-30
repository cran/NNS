#' VN Correlation
#'
#' Returns the nonlinear correlation coefficient based on partial moment quadrants measured by frequency or area.  Degree = 0 is frequency, degree = 1 is area.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param order Controls the level of quadrant partitioning.  Defualts to NULL, but lower levels should be called for large (n).
#' @param degree Defaults to 0 for smaller number of observations
#' @keywords correlation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' VN.cor(x,y)
#' @export

VN.cor = function( x, y, order = NULL,
                   degree= ifelse(length(x)<100,0,1)){




  clpm = numeric(0)
  cupm = numeric(0)
  dlpm = numeric(0)
  dupm = numeric(0)

  if(is.null(order)){
    for (i in 1:floor(log(length(x),4))){
      if(min(nchar(partition.map(x,y,i)$master_part))==i){
        order=i-1}
    } }

  partitioned_df=partition.map(x,y,order)

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


    }}


  nonlin_cor = (sum(clpm) +sum(cupm) -sum(dlpm) -sum(dupm))/(sum(clpm)+sum(cupm)+sum(dlpm)+sum(dupm))

  return(nonlin_cor)


}

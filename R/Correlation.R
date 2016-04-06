#' VN Correlation
#'
#' Returns the nonlinear correlation coefficient based on partial moment quadrants measured by frequency or area.  Degree = 0 is frequency, degree = 1 is area.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param degree Defaults to 0 for smaller number of observations
#' @param order Number of partial moment quadrants to be generated
#' @keywords correlation
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.cor(x,y)}
#' @export

VN.cor = function( x, y,order=ceiling(log10(length(x))),
                   degree= ifelse(length(x)<100,0,1)){

  partitioned_df = partition.map(x, y,order,degree)

  clpm = numeric(0)
  cupm = numeric(0)
  dlpm = numeric(0)
  dupm = numeric(0)

  for(item in unique(partitioned_df$master_part)){
    sub_x = partitioned_df[partitioned_df$master_part == item, 'x']
    sub_y = partitioned_df[partitioned_df$master_part == item, 'y']
    clpm = c(clpm, Co_LPM(degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
    cupm = c(cupm, Co_UPM(degree,mean(sub_x),mean(sub_y), sub_x, sub_y))
    dlpm = c(dlpm, D_LPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
    dupm = c(dupm, D_UPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))


  }

  nonlin_cor = (sum(clpm) +sum(cupm) -sum(dlpm) -sum(dupm))/(sum(clpm)+sum(cupm)+sum(dlpm)+sum(dupm))

  return(nonlin_cor)


}

#' VN Dependence
#'
#' Returns the dependence between two variables based on higher order partial moment correlations measured by frequency or area.
#'
#' @param x Variable 1
#' @param y Variable 2
#' @param degree Defaults to 0 for smaller number of observations
#' @param order Number of partial moment quadrants to be generated
#' @keywords dependence
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.dep(x,y)}
#' @export

VN.dep = function( x, y,order=ceiling(log10(length(x))),
                   degree= ifelse(length(x)<100,0,1)){

  if(order==1){return("Please Increase the Order Specification")}

  partitioned_df = partition.map(x, y,order,degree)

  clpm = numeric(0)
  cupm = numeric(0)
  dlpm = numeric(0)
  dupm = numeric(0)
  rhos = numeric(0)

  for(item in unique(partitioned_df$master_part)){
    sub_x = partitioned_df[partitioned_df$master_part == item, 'x']
    sub_y = partitioned_df[partitioned_df$master_part == item, 'y']
    clpm = c(clpm, Co_LPM(degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
    cupm = c(cupm, Co_UPM(degree,mean(sub_x),mean(sub_y), sub_x, sub_y))
    dlpm = c(dlpm, D_LPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))
    dupm = c(dupm, D_UPM(degree,degree, mean(sub_x),mean(sub_y),sub_x, sub_y))



  }


  for(i in 1:order){

  rhos[i] =  abs((clpm[i]+cupm[i]-dlpm[i]-dupm[i]) / (clpm[i]+cupm[i]+dlpm[i]+dupm[i]))
  }

  plot(x,y)

  m<- rbind(VN.cor(x, y,order,degree),sum(na.omit(rhos))/length(na.omit(rhos)))

  rownames(m) = c("Correlation","Dependence")

  print(m)

  ### Regression Dependence

  return(sum(na.omit(rhos))/length(na.omit(rhos)))

}

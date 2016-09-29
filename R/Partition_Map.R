#' Partition Map
#'
#'  Creates partitions based on quadrant means, assigning observations to those quadrants.  Needed for correlation, dependence, regression routines.  Default degree = 1 for area, but routines have their own conditional degree specifications built in
#' @param x Variable 1
#' @param y Variable 2
#' @param type Controls the partitioning basis.  Set to \code{type="XONLY"} for X-axis based partitioning.  Defaults to NULL for both X and Y-axis partitioning.
#' @param order Number of partial moment quadrants to be generated.
#' @param override Reduces minimum number of necessary observations in a quadrant to 1 when \code{override=TRUE}.
#' @return Returns both a dataframe ("df") of X and Y observations with their partition assignment in the 3d column; and the regression points ("regression.points") for that given \code{order}.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' partition.map(x,y)
#' ## Dataframe of observations and partitions
#' partition.map(x,y,order=1)$df
#' ## Regression points
#' partition.map(x,y,order=1)$regression.points
#' @export

partition.map = function(x, y,type=NULL,order= NULL,override=FALSE){

  temp_df = data.frame(x=x, y=y)
  temp_df[,'temp_part'] = 'p'
  temp_df[,'master_part'] = 'p'

  regression.points = data.frame(matrix(ncol = 2))


  if(override==FALSE){min.obs=4}else{min.obs=1}

  if(is.null(order)){max.order = ceiling(log(length(x),2))}
  else{max.order=order}

### X and Y partition
  if(is.null(type)){
  for(i in 1:max.order){
    regression.points = data.frame(matrix(ncol = 2))

    for(item in unique(temp_df$master_part)){

       if(nchar(item)==i && length(temp_df[temp_df$master_part == item,'x'])>=min.obs){

        tmp_xbar = mean(temp_df[temp_df$master_part == item,'x'])
        tmp_ybar = mean(temp_df[temp_df$master_part == item, 'y'])


        temp_df[temp_df$x > tmp_xbar & temp_df$y > tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x > tmp_xbar & temp_df$y > tmp_ybar & temp_df$master_part == item,'master_part'], 1, sep = '')
        temp_df[temp_df$x <= tmp_xbar & temp_df$y > tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x <= tmp_xbar & temp_df$y > tmp_ybar & temp_df$master_part == item,'master_part'], 2, sep = '')
        temp_df[temp_df$x > tmp_xbar & temp_df$y <= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x > tmp_xbar & temp_df$y <= tmp_ybar & temp_df$master_part == item,'master_part'], 3, sep = '')
        temp_df[temp_df$x <= tmp_xbar & temp_df$y <= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x <= tmp_xbar & temp_df$y <= tmp_ybar & temp_df$master_part == item,'master_part'], 4, sep = '')


            regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

        }
    }

      temp_df[,'master_part'] = temp_df[, 'temp_part']

      if(min(nchar(temp_df[,'master_part'] ))<max(nchar(temp_df[,'master_part'] ))){break}

    }
    q=length(regression.points[,1])
  colnames(regression.points)=c('X','Y')

  return(list("df"=temp_df[, c('x', 'y', 'master_part')],"regression.points"=regression.points[order(regression.points[,1]),][-q,]))
  }

### X only partition
  if(!is.null(type)&&type!="YONLY"){
    for(i in 1:max.order){
      regression.points = data.frame(matrix(ncol = 2))

      for(item in unique(temp_df$master_part)){

        if(nchar(item)==i && length(temp_df[temp_df$master_part == item,'x'])>=min.obs){

          tmp_xbar = mean(temp_df[temp_df$master_part == item,'x'])
          tmp_ybar = mean(temp_df[temp_df$master_part == item, 'y'])


          temp_df[temp_df$x > tmp_xbar  & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x > tmp_xbar  & temp_df$master_part == item,'master_part'], 1, sep = '')
          temp_df[temp_df$x <= tmp_xbar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x <= tmp_xbar  & temp_df$master_part == item,'master_part'], 2, sep = '')


          regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

        }
      }

      temp_df[,'master_part'] = temp_df[, 'temp_part']

      if(min(nchar(temp_df[,'master_part'] ))<max(nchar(temp_df[,'master_part'] ))){break}

    }

    q=length(regression.points[,1])

    colnames(regression.points)=c('X','Y')
    return(list("df"=temp_df[, c('x', 'y', 'master_part')],"regression.points"=regression.points[order(regression.points[,1]),][-q,]))
  }





### Y only partition
if(!is.null(type)&&type=="YONLY"){
  for(i in 1:max.order){
    regression.points = data.frame(matrix(ncol = 2))

    for(item in unique(temp_df$master_part)){

      if(nchar(item)==i && length(temp_df[temp_df$master_part == item,'y'])>=min.obs){

        tmp_xbar = mean(temp_df[temp_df$master_part == item,'x'])
        tmp_ybar = mean(temp_df[temp_df$master_part == item, 'y'])


        temp_df[temp_df$y > tmp_ybar  & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$y > tmp_ybar  & temp_df$master_part == item,'master_part'], 3, sep = '')
        temp_df[temp_df$y <= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$y <= tmp_ybar  & temp_df$master_part == item,'master_part'], 4, sep = '')


        regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

      }
    }

    temp_df[,'master_part'] = temp_df[, 'temp_part']

    if(min(nchar(temp_df[,'master_part'] ))<max(nchar(temp_df[,'master_part'] ))){break}

  }

  q=length(regression.points[,1])

  colnames(regression.points)=c('X','Y')
  return(list("df"=temp_df[, c('x', 'y', 'master_part')],"regression.points"=regression.points[order(regression.points[,2]),][-q,]))
}


}


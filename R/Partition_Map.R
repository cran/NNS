#' NNS Partition Map
#'
#'  Creates partitions based on quadrant means, assigning observations to those quadrants.  Needed for correlation, dependence, regression routines.  Default degree = 1 for area, but routines have their own conditional degree specifications built in
#' @param x Variable 1
#' @param y Variable 2
#' @param Voronoi Displays a Voronoi type diagram using partial moment quadrants.  Defaults to FALSE.
#' @param type Controls the partitioning basis.  Set to \code{type="XONLY"} for X-axis based partitioning.  Defaults to NULL for both X and Y-axis partitioning.
#' @param order Number of partial moment quadrants to be generated.  \code{order="max"} will institute a perfect fit.
#' @param overfit Reduces minimum number of necessary observations in a quadrant to 1 when \code{overfit=TRUE}.  In the instances where \code{"regression.points"} fail to be generated in the output, re-run partitioning with \code{overfit=TRUE} for the given \code{order}.
#' @param noise.reduction \code{noise.reduction="median"} uses medians instead of means for partitions, while \code{noise.reduction="mode"} uses modes instead of means for partitions.  Defaults to \code{noise.reduction="mean"}, while \code{noise.reduction=NULL} will partition quadrant to a single observation for a given \code{order}.
#' @return Returns both a dataframe \code{"df"} of X and Y observations with their partition assignment in the 3d column; and the regression points \code{"regression.points"} for that given \code{order}.
#' @keywords partitioning, cluster
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' NNS.part(x,y)
#' ## Dataframe of observations and partitions
#' NNS.part(x,y,order=1)$df
#' ## Regression points
#' NNS.part(x,y,order=1)$regression.points
#' ## Voronoi style plot
#' NNS.part(x,y,Voronoi=TRUE)
#' @export

NNS.part = function(x, y,Voronoi=FALSE,type=NULL,order= NULL,overfit=FALSE,noise.reduction="mean"){

  if(!is.null(order)){if(order==0) order=1}

  temp_df = data.frame(x=x, y=y)
  temp_df[,'quadrant'] = 'q'
  temp_df[,'obs'] = c(1:length(x))

  mode=function(x) {
    if(length(x)>1){
      d <- density(x)
      d$x[which.max(d$y)]
    }else{x}
  }

  regression.points = data.frame(matrix(ncol = 2))

  if(Voronoi==T){
    plot(x,y,col='steelblue',cex.lab=2,xlab = "X",ylab="Y")}



  if(is.null(order)){
    max.order = ceiling(log(length(y),4))
    overfit=overfit}
  else{
    max.order=order
    if(order==ceiling(log2(length(y)))){
        overfit=T
        type="XONLY"}else{
        overfit=overfit
        type=type
        }}

  if(noise.reduction=='off'){overfit=T}else{overfit=overfit}

  if(overfit==T){min.obs=1}else{min.obs=4}

  ### X and Y partition
  if(is.null(type)){
    for(i in 1:max.order){
      regression.points = data.frame(matrix(ncol = 2))

      for(item in unique(temp_df$quadrant)){
        if(nchar(item)==i){
          sub_df=temp_df[temp_df$quadrant == item,]
          sub_x=sub_df$x
          } else {break}

        if(length(sub_x)>=min.obs){

          sub_y=sub_df$y

          if(noise.reduction=='off'){
            tmp_xbar = mean(sub_x)
            tmp_ybar = mean(sub_y)} else {
              if(noise.reduction=='mode'){
                tmp_xbar = mode(sub_x)
                tmp_ybar = mode(sub_y)} else {
                  if(noise.reduction=='median'){
                    tmp_xbar = median(sub_x)
                    tmp_ybar = median(sub_y)}else{
                      tmp_xbar = mean(sub_x)
                      tmp_ybar = mean(sub_y)}
                }}

          regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

          if(Voronoi==T){
            if(!is.null(noise.reduction)){
              if(noise.reduction=='mean'){
                segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3)
                segments(mean(sub_x),min(sub_y),mean(sub_x),max(sub_y),lty=3)}else{  if(noise.reduction=='mode'){
                  segments(min(sub_x),mode(sub_y),max(sub_x),mode(sub_y),lty=3)
                  segments(mode(sub_x),min(sub_y),mode(sub_x),max(sub_y),lty=3) } else {
                    segments(min(sub_x),median(sub_y),max(sub_x),median(sub_y),lty=3)
                    segments(median(sub_x),min(sub_y),median(sub_x),max(sub_y),lty=3)
                  }
                }
            }
            if(is.null(noise.reduction)){
              segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3)
              segments(mean(sub_x),min(sub_y),mean(sub_x),max(sub_y),lty=3)
            }
          }


          quadrant.1=sub_df[(sub_x > tmp_xbar & sub_y > tmp_ybar),'obs']
          temp_df$quadrant[quadrant.1]=paste(temp_df$quadrant[quadrant.1],1,sep = '')

              sub_df=sub_df[!(sub_df$obs%in%quadrant.1),]
              sub_x=sub_df$x;sub_y=sub_df$y

          quadrant.2=sub_df[(sub_x <= tmp_xbar & sub_y > tmp_ybar),'obs']
          temp_df$quadrant[quadrant.2]=paste(temp_df$quadrant[quadrant.2],2,sep = '')

              sub_df=sub_df[!(sub_df$obs%in%quadrant.2),]
              sub_x=sub_df$x;sub_y=sub_df$y

          quadrant.3=sub_df[(sub_x > tmp_xbar & sub_y <= tmp_ybar),'obs']
          temp_df$quadrant[quadrant.3]=paste(temp_df$quadrant[quadrant.3],3,sep = '')

          quadrant.4=sub_df[!(sub_df$obs%in%quadrant.3),'obs']
          temp_df$quadrant[quadrant.4]=paste(temp_df$quadrant[quadrant.4],4,sep = '')

        }

      }


      if(min(nchar(temp_df[,'quadrant'] ))<max(nchar(temp_df[,'quadrant'] ))){break}

    }

    colnames(regression.points)=c('X','Y')
    q=length(regression.points[,1])
    regression.points=(regression.points[order(regression.points[,1]),][-q,])

    if(Voronoi==T){
      points(regression.points[,1],regression.points[,2],pch=15,lwd=2,col='red')
      title(main=paste(paste0("NNS Order = ",ifelse(is.null(order),i,max.order)),sep="\n"),cex.main=2)
    }

    return(list("df"=temp_df[, c('x', 'y', 'quadrant')],"regression.points"=regression.points))
  }

  ### X only partition
  if(!is.null(type)){
    for(i in 1:max.order){
      regression.points = data.frame(matrix(ncol = 2))

      for(item in unique(temp_df$quadrant)){
        if(nchar(item)==i){
          sub_df=temp_df[temp_df$quadrant == item,]
          sub_x=sub_df$x} else {break}
        if(length(sub_x)>=min.obs){
          sub_y=sub_df$y
          if(noise.reduction=='off'){
            tmp_xbar = mean(sub_x)
            tmp_ybar = mean(sub_y)} else {
              if(noise.reduction=='mode'){
                tmp_xbar = mode(sub_x)
                tmp_ybar = mode(sub_y)} else {
                  if(noise.reduction=='median'){
                    tmp_xbar = median(sub_x)
                    tmp_ybar = median(sub_y)}else{
                      tmp_xbar = mean(sub_x)
                      tmp_ybar = mean(sub_y)}
                }}

          regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

          section.1=sub_df[sub_x > tmp_xbar,'obs']
          temp_df$quadrant[section.1]=paste(temp_df$quadrant[section.1],1,sep = '')

          section.2=sub_df[!(sub_df$obs%in%section.1),'obs']
          temp_df$quadrant[section.2]=paste(temp_df$quadrant[section.2],2,sep = '')

        }


      }



      if(min(nchar(temp_df[,'quadrant'] ))<max(nchar(temp_df[,'quadrant'] ))){break}

    }


    colnames(regression.points)=c('X','Y')
    q=length(regression.points[,1])
    regression.points=(regression.points[order(regression.points[,1]),][-q,])

    if(Voronoi==T){
      points(regression.points[,1],regression.points[,2],pch=15,lwd=2,col='red')
      title(main=paste(paste0("NNS Order = ",ifelse(is.null(order),i,max.order)),sep="\n"),cex.main=2)
    }

    return(list("df"=temp_df[, c('x', 'y', 'quadrant')],"regression.points"=regression.points))
  }


}

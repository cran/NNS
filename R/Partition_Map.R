#' Partition Map
#'
#'  Creates partitions based on quadrant means, assigning observations to those quadrants.  Needed for correlation, dependence, regression routines.  Default degree = 1 for area, but routines have their own conditional degree specifications built in
#' @param x Variable 1
#' @param y Variable 2
#' @param Voronoi Displays a Voronoi type diagram using partial moment quadrants.  Defaults to FALSE.
#' @param type Controls the partitioning basis.  Set to \code{type="XONLY"} for X-axis based partitioning.  Defaults to NULL for both X and Y-axis partitioning.
#' @param order Number of partial moment quadrants to be generated.  \code{order="max"} will institute a perfect fit.
#' @param overfit Reduces minimum number of necessary observations in a quadrant to 1 when \code{overfit=TRUE}.
#' @param noise.reduction \code{noise.reduction="median"} uses medians instead of means for partitions, while \code{noise.reduction="mode"} uses modes instead of means for partitions.  Defaults to \code{noise.reduction="mean"}, while \code{noise.reduction=NULL} will partition quadrant to a single observation for a given \code{order}.
#' @return Returns both a dataframe \code{"df"} of X and Y observations with their partition assignment in the 3d column; and the regression points \code{"regression.points"} for that given \code{order}.
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
#' ## Voronoi style plot
#' partition.map(x,y,Voronoi=TRUE)
#' @export

partition.map = function(x, y,Voronoi=FALSE,type=NULL,order= NULL,overfit=FALSE,noise.reduction="mean"){

  if(!is.null(order)){if(order==0) order=1}

  temp_df = data.frame(x=x, y=y)
  temp_df[,'temp_part'] = 'q'
  temp_df[,'quadrant'] = 'q'

  #temp_df=data.table(temp_df)
  #setkey(temp_df,quadrant)


  mode=function(x) {
    if(length(x)>1){
    d <- density(x)
    d$x[which.max(d$y)]
  }else{x}
  }

  regression.points = data.frame(matrix(ncol = 2))

  if(Voronoi==T){
    plot(x,y,col='steelblue',cex.lab=2,xlab = "X",ylab="Y")}


  if(!is.null(order) | overfit==T){min.obs=1}else{min.obs=4}
  if(is.null(noise.reduction)){min.obs=1}
  if(is.null(order)){max.order = ceiling(log(length(x),2))}
  else{max.order=order}

  ### X and Y partition
  if(is.null(type)){
    for(i in 1:max.order){
      regression.points = data.frame(matrix(ncol = 2))



      for(item in unique(temp_df$quadrant)){
        sub_x=temp_df[temp_df$quadrant == item,'x']
        #sub_x= temp_df[quadrant==item,x]
        if(nchar(item)==i && length(sub_x)>=min.obs){

          sub_y=temp_df[temp_df$quadrant == item,'y']
          #sub_y= temp_df[quadrant==item,y]

          if(is.null(noise.reduction)){
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
            if(is.null(noise.reduction)){
              segments(min(sub_x),mean(sub_y),max(sub_x),mean(sub_y),lty=3)
              segments(mean(sub_x),min(sub_y),mean(sub_x),max(sub_y),lty=3)}else{  if(noise.reduction=='mode'){
                segments(min(sub_x),mode(sub_y),max(sub_x),mode(sub_y),lty=3)
                segments(mode(sub_x),min(sub_y),mode(sub_x),max(sub_y),lty=3) } else {
                  segments(min(sub_x),median(sub_y),max(sub_x),median(sub_y),lty=3)
                  segments(median(sub_x),min(sub_y),median(sub_x),max(sub_y),lty=3)

                }

                }
          }

          quadrant.1=which(temp_df$x > tmp_xbar & temp_df$y > tmp_ybar & temp_df$quadrant == item)
          temp_df$temp_part[quadrant.1]=paste(temp_df$quadrant[quadrant.1],1,sep = '')
        # temp_df[x>tmp_xbar & y>tmp_ybar & quadrant==item, temp_part:= .(temp_part = paste0(quadrant,1)) ]

          quadrant.2=which(temp_df$x <= tmp_xbar & temp_df$y > tmp_ybar & temp_df$quadrant == item)
          temp_df$temp_part[quadrant.2]=paste(temp_df$quadrant[quadrant.2],2,sep = '')
         # temp_df[x<=tmp_xbar & y>tmp_ybar & quadrant==item, temp_part:= .(temp_part= paste0(quadrant,2))]

          quadrant.3=which(temp_df$x > tmp_xbar & temp_df$y <= tmp_ybar & temp_df$quadrant == item)
          temp_df$temp_part[quadrant.3]=paste(temp_df$quadrant[quadrant.3],3,sep = '')
         # temp_df[x>tmp_xbar & y<=tmp_ybar & quadrant==item, temp_part:= .(temp_part = paste0(quadrant,3))]

          quadrant.4=which(temp_df$x <= tmp_xbar & temp_df$y <= tmp_ybar & temp_df$quadrant == item)
          temp_df$temp_part[quadrant.4]=paste(temp_df$quadrant[quadrant.4],4,sep = '')
         # temp_df[x<=tmp_xbar & y<=tmp_ybar & quadrant==item, temp_part:= .(temp_part = paste0(quadrant,4))]


        }

      }

      #temp_df[, quadrant:= .(quadrant=temp_part)]
      temp_df$quadrant=temp_df$temp_part
      if(min(nchar(temp_df[,4] ))<max(nchar(temp_df[,4] ))){break}

    }
    q=length(regression.points[,1])
    colnames(regression.points)=c('X','Y')
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
        sub_x=temp_df[temp_df$quadrant == item,'x']
        if(nchar(item)==i && length(sub_x)>=min.obs){
          sub_y=temp_df[temp_df$quadrant == item,'y']
          if(is.null(noise.reduction)){
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

          section.1=which(temp_df$x > tmp_xbar  & temp_df$quadrant == item)
          temp_df$temp_part[section.1]=paste(temp_df$quadrant[section.1],1,sep = '')

          section.2=which(temp_df$x <= tmp_xbar  & temp_df$quadrant == item)
          temp_df$temp_part[section.2]=paste(temp_df$quadrant[section.2],2,sep = '')

        }


      }

      temp_df[,'quadrant'] = temp_df[, 'temp_part']

      if(min(nchar(temp_df[,'quadrant'] ))<max(nchar(temp_df[,'quadrant'] ))){break}

    }

    q=length(regression.points[,1])

    colnames(regression.points)=c('X','Y')
    regression.points=(regression.points[order(regression.points[,1]),][-q,])
    return(list("df"=temp_df[, c('x', 'y', 'quadrant')],"regression.points"=regression.points))
  }


}


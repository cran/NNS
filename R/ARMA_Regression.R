#' VN.ARMA.reg
#'
#' VN.reg function with reduced output specifically for VN.ARMA routine
#'
#'
#' @param x Independent Variable
#' @param y Dependent Variable
#' @param order Number of partial moment quadrants to be generated
#' @param point.est Value to be fitted
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.ARMA.reg(x,y)}


VN.ARMA.reg = function (x, y,
                          order=max(2,ceiling(log10(length(x)))),
                          point.est = NULL){

  temp_df = data.frame(x=x, y=y)
  temp_df[,'temp_part'] = 'p'
  temp_df[,'master_part'] = 'p'

  regression.points = data.frame(matrix(ncol = 2))
  Regression.Coefficients = data.frame(matrix(ncol=3))

  names(Regression.Coefficients) = c('Coefficient','X Lower Range','X Upper Range')



  if(order==1){return("Please Increase the Order Specification")}

  if(order >1){
    for(i in 1:(order-1)){

      for(item in unique(temp_df$master_part)){
        tmp_xbar = mean(temp_df[temp_df$master_part == item,'x'])
        tmp_ybar = mean(temp_df[temp_df$master_part == item, 'y'])



        temp_df[temp_df$x >= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'master_part'], 1, sep = '')
        temp_df[temp_df$x <= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x <= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'master_part'], 2, sep = '')
        temp_df[temp_df$x >= tmp_xbar & temp_df$y <= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar & temp_df$y <= tmp_ybar & temp_df$master_part == item,'master_part'], 3, sep = '')
        temp_df[temp_df$x <= tmp_xbar & temp_df$y <= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x <= tmp_xbar & temp_df$y <= tmp_ybar & temp_df$master_part == item,'master_part'], 4, sep = '')

        if(nchar(item)==order-1){

          regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

        }


      }

      temp_df[,'master_part'] = temp_df[, 'temp_part']
    }


  }



  ###Endpoints
  x0 = temp_df[order(temp_df$x),][1,2]
  x.max = temp_df[order(temp_df$x),][length(x),2]

  regression.points[1,2] = x0
  regression.points[1,1] = min(x)

  regression.points[length(regression.points[,2])+1,2] = x.max
  regression.points[length(regression.points[,1]),1] = max(x)



  ###Regression Equation

  regression.points = na.omit(regression.points[order(regression.points),])


  q=length(regression.points[,1])

  for(i in 1:q){

    rise = regression.points[i+1,2] - regression.points[i,2]
    run = regression.points[i+1,1] - regression.points[i,1]

    Regression.Coefficients[i,] = cbind((rise/run),regression.points[i,1],regression.points[i+1,1])
    Regression.Coefficients[q,] = cbind(1,regression.points[i,1],regression.points[i,1]+1e-15)
  }

  Regression.Coefficients= na.omit(Regression.Coefficients)

  ### Fitted Values
  p = length((Regression.Coefficients)[,1])


  fitted = numeric()
  fitted.new = numeric()

  for (i in 1:p){

    z=(which(x>=Regression.Coefficients[i,2] & x<Regression.Coefficients[(i),3]))

    z.diff = ((x[z]- Regression.Coefficients[i,2])*Regression.Coefficients[i,1])+regression.points[i,2]

    if(is.null(point.est)){point.est.y = NULL} else{

      if(!is.null(point.est) && point.est>=Regression.Coefficients[i,2] && point.est<Regression.Coefficients[i,3]){ point.est.y = (point.est - Regression.Coefficients[i,2])*(Regression.Coefficients[i,1])+regression.points[i,2]}

      else{if(!is.null(point.est) && point.est<Regression.Coefficients[1,2]){
        point.est.y = ((point.est - Regression.Coefficients[2,2])*(Regression.Coefficients[1,1]))+(regression.points[2,2])
      }

        else{if(!is.null(point.est) && point.est>Regression.Coefficients[p,2]){point.est.y = ((point.est - Regression.Coefficients[(p-1),2])*(Regression.Coefficients[(p-1),1]))+(regression.points[(p-1),2])
        }
        }
      }
    }


    fitted.new =  cbind(z,z.diff)


    fitted = rbind(fitted,fitted.new)
    fitted = fitted[order(fitted[,1]),]

  }




  if(!is.null(point.est)) {return(point.est.y)}

}

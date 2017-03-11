#' NNS ARMA
#'
#' Autoregressive model incorporating nonlinear regressions of component series.
#'
#' @param variable a numeric vector.
#' @param h integer; 1 (default) Number of periods to forecast.
#' @param Training.set \code{NULL} (defualt) numeric; Sets the number of variable observations \code{(variable[1:Training.set])} to monitor performance of forecast over in-sample range.
#' @param Seasonal.Factor logical or integer; \code{TRUE} (default) Automatically selects the best seasonal lag from the seasonality test.  To use weighted average of all seasonal lags set to \code{(Seasonal.Factor=FALSE)}.  Otherwise, directly input known frequency integer lag to use, i.e. \code{(Seasonal.Factor=12)} for monthly data.
#' @param Negative.Values logical; \code{FALSE} (default) If the variable can be negative, set to \code{(Negative.Values=TRUE)}.
#' @param Method options:("lin","nonlin","both"); To select the regression type of the component series, defaults to \code{(Method="both")} where both linear and nonlinear estimates are generated.  To use a nonlineaer regression, set to \code{(Method="nonlin")}; to use a linear regression set to \code{Method="lin"}.
#' @param Dynamic logical; \code{FALSE} (default) To update the seasonal factor with each forecast point, set to \code{(Dynamic=TRUE)}.  The default is \code{(Dynamic=FALSE)} to retain the original seasonal factor from the inputted variable for all ensuing \code{h}.
#' @param stn numeric [0,1]; Signal to noise parameter, sets the threshold of \code{(NNS.dep)} which reduces \code{("order")} when \code{(order=NULL)}.  Defaults to \code{(stn=0)} to allow for maximum fit.
#' @param plot logical; \code{TRUE} (default) Returns the plot of all periods exhibiting seasonality and the variable level reference in upper panel.  Lower panel returns original data and forecast.
#' @param Seasonal_plot logical; \code{TRUE} (default) Adds the seasonality plot above the forecast.  Will be set to \code{FALSE} if no seasonality is detected or \code{Seasonal.Factor} is set to an integer value.
#' @param Intervals logical; \code{FALSE} (default) Plots the surrounding forecasts around the final estimate when \code{(Intervals=TRUE)} and \code{(Seasonal.Factor=FALSE)}.  There are no other forecasts to plot when a single \code{Seasonal.Factor} is selected.
#' @return Returns a vector of forecasts of length \code{(h)}.
#' @note \code{(Seasonal.Factor=FALSE)} can be a very comutationally expensive exercise due to the number of seasonal periods detected.
#' @keywords Autoregressive model
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' ## Nonlinear NNS.ARMA using AirPassengers monthly data and 12 period lag
#' NNS.ARMA(AirPassengers,h=45,Training.set=100,Seasonal.Factor=12,Method='nonlin')
#' @export



# Autoregressive Model
NNS.ARMA <- function(variable,h=1,Training.set = NULL, Seasonal.Factor = TRUE ,Negative.Values = FALSE, Method = "both", Dynamic = FALSE,stn=0,plot=TRUE,Seasonal_plot=TRUE,Intervals=FALSE){

  if(Intervals==TRUE && is.numeric(Seasonal.Factor)){stop('Hmmm...Seems you have "Intervals" and "Seasonal.Factor" selected.  Please set "Intervals=F" or "Seasonal.Factor=F"')}

  if(Intervals==TRUE && Seasonal.Factor==TRUE){stop('Hmmm...Seems you have "Intervals" and "Seasonal.Factor" selected.  Please set "Intervals=F" or "Seasonal.Factor=F"')}


  variable=as.numeric(variable)
  OV = variable

  if(!is.null(Training.set)){
    variable = variable[1:Training.set]
    FV = variable[1:Training.set]
  } else {
    Training.set = length(variable)
    variable = variable
    FV = variable
  }

  Estimates = numeric()

  # Weight and lag function for seasonality...
  ARMA.seas.weighting=function(){
    if(is.null(ncol(M))){
      return(list(lag=M[1],Weights=1))}

    if(ncol(M)==1){
      return(list(lag=1,Weights=1))}

    if(ncol(M)>1){
        if(Seasonal.Factor==TRUE){
            lag = seas.matrix$best.period
            Weights=1
            return(list(lag=lag,Weights=Weights))
        }

      Observation.sum = sum(1/sqrt(M[,1]))
      Observation.weighting = (1/sqrt(M[,1]))

      Lag.sum = sum(M[,3]-M[,2])
      Lag.weighting = (M[,3]-M[,2])


      Weights = (Lag.weighting*Observation.weighting) / sum(Lag.weighting*Observation.weighting)

      # Determine lag from seasonality test
      if(Seasonal.Factor==FALSE){lag<- M[,1]}
      if(is.numeric(Seasonal.Factor)){
        lag<- Seasonal.Factor
        Weights=1}
      return(list(lag=lag,Weights=Weights))
    }}

  #Vectors generator for 1:lag
  generate.vectors=function(lag){
    Component.series = list()
    Component.index = list()

    for (i in 1:length(lag)){
      Component.series[[paste('Series.',i,sep="")]] <- na.omit(rev(variable[seq(aa+1,1,-lag[i])]))
      Component.index[[paste('Index.',i,sep="")]] <- (1:length(na.omit(Component.series[[i]])))
    }
    return(list(Component.index=Component.index,Component.series=Component.series))
  }



  if(Dynamic == TRUE){
    for (j in 0:(h-1)){
      if(is.numeric(Seasonal.Factor)){
        M<-t(Seasonal.Factor)
        lag=Seasonal.Factor
        Weights=1
        Seasonal_plot=F
      } else {

        seas.matrix = NNS.seas(variable,plot=F)
        if(!is.list(seas.matrix)){
          M<- t(seas.matrix)} else {
            M<- seas.matrix$all.periods}
        ASW=ARMA.seas.weighting()
        lag=ASW$lag
        Weights=ASW$Weights
      }

      a=length(FV)
      aa=length(variable)

      # Generate vectors for 1:lag
      GV=generate.vectors(lag)
      Component.index=GV$Component.index
      Component.series=GV$Component.series
      # Regression on Component Series
      Regression.Estimates = numeric()
      Coefficients = numeric()
      Estimate.band= list()


      if(Method=='nonlin' | Method=='both'){
        if(Negative.Values==TRUE){
          Regression.Estimates=sapply(1:length(lag),function(i) NNS.reg(Component.index[[i]],Component.series[[i]],point.est = (length(Component.series[[i]])+1),return.values = TRUE,order = NULL,plot = FALSE,stn=stn)$Point.est)
        }
        if(Negative.Values==FALSE){
          Regression.Estimates=sapply(1:length(lag),function(i) NNS.reg(Component.index[[i]],Component.series[[i]],point.est = (length(Component.series[[i]])+1),return.values = TRUE,order = NULL,plot = FALSE,stn=stn)$Point.est)
          Regression.Estimates=pmax(0,Regression.Estimates)
        }
        Nonlin.estimates=sum(Regression.Estimates*Weights)
      }#Linear == F

      if(Method=='lin' | Method=='both'){
        if(Negative.Values==FALSE){
          Regression.Estimates=sapply(1:length(lag),function(i) coef(lm(Component.series[[i]]~Component.index[[i]]))[1]+(coef(lm(Component.series[[i]]~Component.index[[i]]))[2]*(length(Component.series[[i]])+1)))
          Regression.Estimates=pmax(0,Regression.Estimates)
        }
        if(Negative.Values==TRUE){
          Regression.Estimates=sapply(1:length(lag),function(i) coef(lm(Component.series[[i]]~Component.index[[i]]))[1]+(coef(lm(Component.series[[i]]~Component.index[[i]]))[2]*(length(Component.series[[i]])+1)))
        }
        Lin.estimates=sum(Regression.Estimates*Weights)
      }#Linear==T

      if(Intervals==TRUE){
        if(Method=='both'){
          Estimate.band[[j+1]]=c(Nonlin.estimates,Lin.estimates)}
        if(Method=='nonlin'){
          Estimate.band[[j+1]]=Nonlin.estimates}
        if(Method=='lin'){
          Estimate.band[[j+1]]=Lin.estimates}
      }


      if(Method=='both'){Estimates[j+1]=mean(c(Lin.estimates,Nonlin.estimates))}else{
        Estimates[j+1] = sum(Regression.Estimates*Weights)}
      variable = c(variable,(Estimates[j+1]))
      FV=variable

    } #j loop
  }  #Dynamic {}

  else {
    if(is.numeric(Seasonal.Factor)){
      M<-t(Seasonal.Factor)
      lag=Seasonal.Factor
      Weights=1
      Seasonal_plot=F
    } else {

    seas.matrix = NNS.seas(variable,plot=F)
    if(!is.list(seas.matrix)){
      M<- t(seas.matrix)} else {
        M<- seas.matrix$all.periods}
    ASW=ARMA.seas.weighting()
    lag=ASW$lag
    Weights=ASW$Weights
    }

    Estimate.band= list()

    for (j in 0:(h-1)){
      a=length(FV)
      aa=length(variable)

      # Generate vectors for 1:lag
      GV=generate.vectors(lag)
      Component.index=GV$Component.index
      Component.series=GV$Component.series

      # Regression on Component Series
      Regression.Estimates = numeric()
      Coefficients = numeric()


      if(Method=='nonlin' | Method=='both'){
        if(Negative.Values==TRUE){
          Regression.Estimates=sapply(1:length(lag),function(i) NNS.reg(Component.index[[i]],Component.series[[i]],point.est = (length(Component.series[[i]])+1),return.values = TRUE,order = NULL,plot = FALSE,stn=stn)$Point.est)
          NL.Regression.Estimates=Regression.Estimates
        }

        if(Negative.Values==FALSE){
          Regression.Estimates=sapply(1:length(lag),function(i) NNS.reg(Component.index[[i]],Component.series[[i]],point.est = (length(Component.series[[i]])+1),return.values = TRUE,order = NULL,plot = FALSE,stn=stn)$Point.est)
          Regression.Estimates=pmax(0,Regression.Estimates)
          NL.Regression.Estimates=Regression.Estimates
        }

        Nonlin.estimates=sum(Regression.Estimates*Weights)
      }#Linear == F

      if(Method=='lin' | Method=='both'){
        if(Negative.Values==FALSE){
          Regression.Estimates=sapply(1:length(lag),function(i) coef(lm(Component.series[[i]]~Component.index[[i]]))[1]+(coef(lm(Component.series[[i]]~Component.index[[i]]))[2]*(length(Component.series[[i]])+1)))
          Regression.Estimates=pmax(0,Regression.Estimates)
          L.Regression.Estimates=Regression.Estimates
        }

        if(Negative.Values==TRUE){
          Regression.Estimates=sapply(1:length(lag),function(i) coef(lm(Component.series[[i]]~Component.index[[i]]))[1]+(coef(lm(Component.series[[i]]~Component.index[[i]]))[2]*(length(Component.series[[i]])+1)))
          L.Regression.Estimates=Regression.Estimates
        }

        Lin.estimates=sum(Regression.Estimates*Weights)

      }#Linear==T

      if(Intervals==TRUE){
        if(Method=='both'){
          Estimate.band[[j+1]]=c(NL.Regression.Estimates,L.Regression.Estimates)}
        if(Method=='nonlin'){
          Estimate.band[[j+1]]=NL.Regression.Estimates}
        if(Method=='lin'){
          Estimate.band[[j+1]]=L.Regression.Estimates}
      }

      if(Method=='both'){
        Estimates[j+1]=mean(c(Lin.estimates,Nonlin.estimates))}else{
          Estimates[j+1] = sum(Regression.Estimates*Weights)}
      variable = c(variable,(Estimates[j+1]))
      FV=variable

    }}  # ELSE {}
  #return(Estimate.band)
  #### PLOTTING
  if(plot==T){
    if(Seasonal_plot==T){
    if(ncol(M)>1){
      par(mfrow=c(2,1))
      plot(M[,Period],M[,Coefficient.of.Variance],
           xlab="Period", ylab="Coefficient of Variance", main = "Seasonality Test",         ylim = c(0,2*abs(sd(FV)/mean(FV))))

      points(M[1,Period],M[1,Coefficient.of.Variance],pch=19,col='red')

      abline(h=abs(sd(FV)/mean(FV)), col="red",lty=5)
      text(mean(M[,Period]),abs(sd(FV)/mean(FV)),pos=3,"Variable Coefficient of Variance",col='red')
    } else {
      par(mfrow=c(2,1))
      plot(1,1,pch=19,col='blue', xlab="Period", ylab="Coefficient of Variance", main = "Seasonality Test",         ylim = c(0,2*abs(sd(FV)/mean(FV))))

      text(1,abs(sd(FV)/mean(FV)),pos=3,"NO SEASONALITY DETECTED",col='red')
    }
}
    plot( OV, type = 'l',lwd=2,main = "Forecast",col='steelblue',  xlim=c(1,max((a+h),length(OV))), ylab="Variable", ylim=c(min(Estimates, OV),max( OV,Estimates)))

    if(Intervals==T){
      for(i in 1:h){
        ys=unlist(Estimate.band[[i]])
        points(rep(Training.set+i,length(ys)),ys,col=rgb(1, 0, 0, 0.0125),pch=15)
      }
      lines((Training.set+1):(Training.set+h),Estimates,type = 'l',lwd=2,lty=1,col='red')

      segments(Training.set,FV[Training.set],Training.set+1,Estimates[1],lwd=2,lty=1,col='red')
      legend('topleft',bty='n',legend = c("Original", paste("Forecast ",h," period(s)",sep = "")),lty = c(1,1),col=c('steelblue','red'),lwd=2)

    } else{

      if(Training.set[1]<length(OV)){
        lines((Training.set+1):(Training.set+h),Estimates,type = 'l',lwd=2,lty=3,col='red')

        segments(Training.set,FV[Training.set],Training.set+1,Estimates[1],lwd=2,lty=3,col='red')
        legend('topleft',bty='n',legend = c("Original", paste("Forecast ",h," period(s)",sep = "")),lty = c(1,2),col=c('steelblue','red'),lwd=2)

      } else {
        lines((Training.set+1):(Training.set+h),Estimates,type = 'l',lwd=2,lty=1,col='red')

        segments(Training.set,FV[Training.set],Training.set+1,Estimates[1],lwd=2,lty=1,col='red')
        legend('topleft',bty='n',legend = c("Original", paste("Forecast ",h," period(s)",sep = "")),lty = c(1,1),col=c('steelblue','red'),lwd=2)
      }


    }
    points(Training.set,variable[Training.set],col="green",pch=18)
    points(Training.set+h,sum(Regression.Estimates*Weights),col="green",pch=18)

    par(mfrow=c(1,1))
  }
  return(Estimates)

}

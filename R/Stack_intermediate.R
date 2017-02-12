#' NNS.stack.intermediate (INTERNAL CALL FOR \link{NNS.stack})
#'
#' Called by \code{NNS.stack} for precision testing in cross-validation.
#'
#' @param IVs.train Training set of independent variables.
#' @param IVs.test Test set of independent variables.
#' @param DV.train Training set of dependent variable.
#' @param DV.test Test set of dependent variable.
#' @param CV.size Sets the cross-validation size if \code{DV.test=NULL}.  Defaults to 0.2 or a 20 percent random sampling of the training set.
#' @param weight Set \code{weight="MSE"} for optimum parameters and weighting based on each base model's \code{"MSE"}.  \code{weight="Feautures"} uses a weighting based on the number of features present, whereby logistic \link{NNS.reg} and \link{NNS.Feature.prob} receive higher weights for more independent variables.  Defaults to \code{"MSE"}.
#' @param text If performing a text classification, set \code{text=TRUE}.  Defaults to FALSE.
#' @param precision Increases speed of computation at the expense of precision.  2 settings offered: \code{"LOW"} ,and \code{"HIGH"}.  \code{"HIGH"} is the limit condition of every observation as a regression point.  \code{precision=NULL} (Defualt) compares both precision types and then returns the best model parameters.
#' @param method Select the NNS method to include in stack.  \code{method=1} selects \link{NNS.reg}; \code{method=2} selects \link{NNS.reg} dimension reduction regression; \code{method=3} selects \link{NNS.Feature.prob}.  Defaults to \code{method=c(1,2,3)}, including all 3 NNS methods in the stack.
#' @param threshold  Sets the correlation threshold for independent variables in \link{NNS.reg}.  Defaults to \code{threshold=0}.
#' @return Returns a vector of fitted values for the dependent variable test set for all models.  \code{"NNS.reg.n.best"} returns the optimum \code{"n.best"} paramater for the \link{NNS.reg} multivariate regression.  \code{"NNS.logistic.order"} returns the optimum \code{"order"} from the \link{NNS.reg} logistic regression.  \code{"reg"} returns \link{NNS.reg} output, \code{"logistic"} returns \link{NNS.reg} logistic regression output, \code{"Feature.prob"} returns \link{NNS.Feature.prob} output, and \code{"stack"} returns the output of the stacked model.  \code{"CV.test"} returns the random cross-validation dependent variable used and \code{"MSE"} returns the MSE from cross-validation.
#' @author Fred Viole, OVVO Financial Systems




NNS.stack.intermediate <- function(IVs.train,DV.train,IVs.test=NULL,DV.test=NULL,CV.size=.2,weight="MSE",text=F,precision=NULL,method=c(1,2,3),threshold = 0){

  IVs.train<- apply(IVs.train,2,as.numeric)
  DV.train<- as.numeric(DV.train)

  n<- ncol(IVs.train)
  n.test<- ncol(IVs.test)

  if(!is.null(IVs.test)){
    if(is.null(n.test)){
    IVs.test= as.numeric(IVs.test)
    } else {

    IVs.test<- apply(IVs.test,2,as.numeric)
    }
    }

  if(!is.null(DV.test)){
    DV.test<- as.numeric(DV.test)
    CV.DV.test<- DV.test}


  l=length(IVs.train[,1])

  # Multivariate Regression Output

  if(is.null(DV.test)){

    test.set=sample(1:l,as.integer(CV.size*l),
                    replace = FALSE)

    CV.IVs.train<- IVs.train[c(-test.set),]
    CV.IVs.test<- IVs.train[c(test.set),]
    CV.DV.train<- DV.train[c(-test.set)]
    CV.DV.test<- DV.train[c(test.set)]

    IVs.train<- CV.IVs.train
    if(is.null(IVs.test)){IVs.test<- CV.IVs.test}
    DV.train<- CV.DV.train
    if(is.null(DV.test)){DV.test<- CV.DV.test}
  }

    if(1%in%method){
     nns.cv=numeric()
      if(is.null(DV.test)){
    for(i in 1:(2*n)){
      nns.cv[i]=mean((NNS.reg(CV.IVs.train, CV.DV.train,precision=precision,
                              point.est = CV.IVs.test, plot=F, n.best = i,text=text,norm='NNS')$Point.est-CV.DV.test)^2)
    }

    nns.1=NNS.reg(IVs.train, DV.train,precision=precision,
                  point.est = IVs.test, plot=F, n.best = which.min(nns.cv),text=text,norm='NNS')$Point.est
      } else {

    for(i in 1:(2*n)){
      nns.cv[i]=mean((NNS.reg(IVs.train, DV.train,precision=precision,
                              point.est = IVs.test, plot=F, n.best = i,text=text,norm='NNS')$Point.est-DV.test)^2)
    }

    nns.1=NNS.reg(IVs.train, DV.train,precision=precision,
                  point.est = IVs.test, plot=F, n.best = which.min(nns.cv),text=text,norm='NNS')$Point.est
      }

    } else {nns.1=0
  nns.cv=1e-10}

  # Logistic Regression Output
  if(2%in%method){
  nns.ord=numeric()
  if(is.null(DV.test)){
    for(i in 1:log(l,2)){
      nns.ord[i]=mean((NNS.reg(CV.IVs.train, CV.DV.train,point.est = CV.IVs.test,plot=F, order = i,text=text, type = "CLASS",dep.order = 1,threshold = threshold)$Point.est-CV.DV.test)^2)
    }

    nns.2=NNS.reg(IVs.train, DV.train,
                  point.est = IVs.test, type = "CLASS",
                  plot=F,order=which.min(nns.ord),dep.order = 1,threshold = threshold)$Point.est

  } else {
    for(i in 1:log(l,2)){
      nns.ord[i]=mean((NNS.reg(IVs.train, DV.train,point.est = IVs.test,
                               plot=F, order = i,text=text, type = "CLASS",dep.order=1,threshold = threshold)$Point.est-DV.test)^2)
    }


    nns.2=NNS.reg(IVs.train, DV.train,
                  point.est = IVs.test, type = "CLASS",
                  plot=F,order=which.min(nns.ord),dep.order=1,threshold = threshold)$Point.est
  }
  } else {nns.2=0
  nns.ord=1e-10}

  # NNS Feature Probability Output
  if(3%in%method){
  if(is.null(DV.test)){
    nns.fp=mean((NNS.Feature.prob(cbind(DV.train,IVs.train), 1,
                                  point.est=CV.IVs.test)$Point.est-CV.DV.test)^2)}         else {nns.fp=mean((NNS.Feature.prob(cbind(DV.train,IVs.train), 1,
                                point.est=IVs.test)$Point.est-DV.test)^2)
                                  }

  nns.3=NNS.Feature.prob(cbind(DV.train,IVs.train), 1,
                         point.est=IVs.test)$Point.est
  } else {nns.3=0
  nns.fp=1e-10}

  ### Weights for combining NNS techniques
  if(weight=="Features"){weights=c(n,n-1,n-2)}
  if(weight=="MSE"){weights=c(max(1e-10,1/min(nns.cv)),max(1e-10,1/min(nns.ord)),max(1e-10,1/nns.fp))}

  weights=pmax(weights,c(0,0,0))
  weights[!(c(1,2,3)%in%method)]<- 0
  weights=weights/sum(weights)



  mse=mean((rowMeans(cbind(weights[1]*nns.1+weights[2]*nns.2+weights[3]*nns.3))-DV.test)^2)

  stack.output=rowMeans(cbind(weights[1]*nns.1+weights[2]*nns.2+weights[3]*nns.3))

  return(list(NNS.reg.n.best=which.min(nns.cv),NNS.logistic.order=which.min(nns.ord),MSE=mse,reg=nns.1,logistic=nns.2,Feature.prob=nns.3,stack=stack.output,CV.test=CV.DV.test))

}

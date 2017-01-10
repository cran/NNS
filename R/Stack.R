#' NNS.stack
#'
#' Prediction model using the predictions of the NNS base models (\link{NNS.reg}, \link{NNS.Feature.prob}) as features (i.e. meta features) for the stacked model.
#'
#' @param IVs.train Training set of independent variables.
#' @param IVs.test Test set of independent variables.
#' @param DV.train Training set of dependent variable.
#' @param DV.test Test set of dependent variable.
#' @param CV.size Sets the cross-validation size if \code{DV.test=NULL}.  Defaults to 0.2 or a 20 percent random sampling of the training set.
#' @param weight Set \code{weight="MSE"} for optimum parameters and weighting based on each base model's \code{"MSE"}.  \code{weight="Feautures"} uses a weighting based on the number of features present, whereby logistic \link{NNS.reg} and \link{NNS.Feature.prob} receive higher weights for more independent variables.  Defaults to \code{"MSE"}.
#' @param text If performing a text classification, set \code{text=TRUE}.  Defaults to FALSE.
#' @param precision Increases speed of computation at the expense of precision.  2 settings offered: \code{"LOW"} ,and \code{"HIGH"}.  \code{"HIGH"} is the limit condition of every observation as a regression point.  \code{precision=NULL} (Defualt) compares both precision types and then returns the best model parameters.
#' @return Returns a vector of fitted values for the dependent variable test set for all models.  \code{"NNS.reg.n.best"} returns the optimum \code{"n.best"} paramater for the \link{NNS.reg} multivariate regression.  \code{"NNS.logistic.order"} returns the optimum \code{"order"} from the \link{NNS.reg} logistic regression.  \code{"reg"} returns \link{NNS.reg} output, \code{"logistic"} returns \link{NNS.reg} logistic regression output, \code{"Feature.prob"} returns \link{NNS.Feature.prob} output, and \code{"stack"} returns the output of the stacked model.  \code{"CV.test"} returns the random cross-validation dependent variable used and \code{"MSE"} returns the MSE from cross-validation.
#' @author Fred Viole, OVVO Financial Systems
#' @note If character variables are used, transform them first to factors using \link{as.factor}, or \link{data.matrix} to ensure overall dataset is numeric.  A multifunction \link{sapply} can also be applied to the overall dataset: \code{data<- sapply(data,function(x){as.factor(x);as.numeric(x)})}.  Then run \code{NNS.stack} with transormed variables.
#' @examples
#'  ## Using 'iris' dataset where predictive attributes are columns 1:4, and the class is column 5.
#'  NNS.stack(iris[,1:4],iris[,5],precision="LOW")
#' @export

NNS.stack <- function(IVs.train,DV.train,IVs.test=NULL,DV.test=NULL,CV.size=.2,weight="MSE",text=F,precision=NULL){

  if(is.null(precision)){
  stack.high.precision = NNS.stack.intermediate(IVs.train=IVs.train,DV.train=DV.train,IVs.test=IVs.test,DV.test=DV.test,CV.size=CV.size,weight=weight,text=text,precision="HIGH")
  stack.low.precision = NNS.stack.intermediate(IVs.train=IVs.train,DV.train=DV.train,IVs.test=IVs.test,DV.test=DV.test,CV.size=CV.size,weight=weight,text=text,precision="LOW")

  low.mse=stack.low.precision$MSE
  high.mse=stack.high.precision$MSE

  if(is.na(stack.low.precision$MSE)){low.mse=0}
  if(is.na(stack.high.precision$MSE)){high.mse=0}

  if(low.mse <= high.mse){
        return(stack.low.precision)
        } else {
          return(stack.high.precision)
        }
  } else {
    NNS.stack.intermediate(IVs.train=IVs.train,DV.train=DV.train,IVs.test=IVs.test,DV.test=DV.test,CV.size=CV.size,weight=weight,text=text,precision=precision)
  }
}

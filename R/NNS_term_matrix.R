#' NNS Term Matrix
#'
#' Generates a term matrix for text classification use in \link{VN.reg}.
#'
#' @param x Text A two column dataset should be used.  Concatenate text from original sources to comply with format.  Also note the possiblity of factors in \code{"DV"}, so \code{"as.numeric(as.character(...))"} is used to avoid issues.
#' @param oos Out-of-sample text dataset to be classified
#' @return Returns the text as independent variables \code{"IV"} and the classification as the dependent variable \code{"DV"}.  Out-of-sample independent variables are returned with \code{"OOS"}.
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' x<- data.frame(cbind(c("sunny","rainy"),c(1,-1)))
#' NNS.term.matrix(x)
#'
#' ### Concatenate Text with space seperator, cbind with "DV"
#' x<- data.frame(cbind(c("sunny","rainy"),c("windy","cloudy"),c(1,-1)))
#' x<- data.frame(cbind(paste(x[,1],x[,2],sep=" "),as.numeric(as.character(x[,3]))))
#' NNS.term.matrix(x)
#' @export


NNS.term.matrix <- function(x, oos=NULL){

  p=length(oos)
  x=t(t(x))
  n=nrow(x)
  unique.vocab=unique(unlist(strsplit(as.character(x[,1]), " ", fixed = TRUE)))

  NNS.TM=data.frame()
  OOS.TM=data.frame()

  for(vocab in unique.vocab){
      for(i in 1:n){
          NNS.TM[i,which(unique.vocab==vocab)]=sum(unlist(strsplit(as.character(x[i,1]), " ", fixed = TRUE))==vocab)
      }
      for(i in 1:p){
          OOS.TM[i,which(unique.vocab==vocab)]=sum(unlist(strsplit(as.character(oos[i]), " ", fixed = TRUE))==vocab)
      }
  }
  NNS.TM=cbind(NNS.TM,x[,2])
  colnames(NNS.TM)=c(unique.vocab,colnames(x)[2])

  return(list("IV"=NNS.TM[,1:length(unique.vocab)],"DV"=as.numeric(as.character(NNS.TM[,length(unique.vocab)+1])),"OOS"=OOS.TM))


}

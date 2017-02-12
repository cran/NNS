#' NNS Dependence Matrix (INTERNAL CALL FOR \link{NNS.dep})
#'
#' Returns the dependence matrix for \link{NNS.dep}.
#'
#' @param x Variable 1 matrix
#' @param order Controls the level of quadrant partitioning.  Default to \code{order=NULL}.  Errors can generally be rectified by setting \code{order=1}.
#' @param degree Defaults to NULL to allow number of observations to be \code{degree} determinant.
#' @return Returns the \code{"Correlation"} and \code{"Dependence"}
#' @keywords dependence, correlation
#' @author Fred Viole, OVVO Financial Systems

NNS.dep.matrix = function(x, order = NULL,
                   degree= NULL){

n= ncol(x)
if(is.null(n)){stop("supply both 'x' and 'y' or a matrix-like 'x'")}

rhos = data.frame()
deps = data.frame()

for(j in 0:(n-1)){
  for(i in 1:(n-j)){

    if((i+j)==i){
      rhos[i+j,i]=1;deps[i+j,i]=1
      rhos[i,i+j]=1;deps[i,i+j]=1
    } else {

      rho=NNS.dep(x[,i],x[,i+j],print.map = F,order=order)$Correlation
      dep=NNS.dep(x[,i],x[,i+j],print.map = F,order=order)$Dependence
      rhos[i+j,i]=rho
      rhos[i,i+j]=rho
      deps[i+j,i]=dep
      deps[i,i+j]=dep
      }
  }
}
colnames(rhos) = colnames(x);colnames(deps) = colnames(x)
rownames(rhos) = colnames(x);rownames(deps) = colnames(x)

return(list("Correlation"=rhos,"Dependence"=deps))
}




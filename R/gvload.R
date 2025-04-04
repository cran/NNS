# Import calls and globalvariable calls

#' @importFrom grDevices adjustcolor rainbow rgb
#' @importFrom graphics abline boxplot legend lines par plot points segments text matplot title axis mtext barplot hist strwidth polygon
#' @importFrom quantmod getSymbols
#' @importFrom Rfast colmeans rowmeans rowsums comb_n
#' @importFrom stats coef cor cov lm na.omit sd median complete.cases resid uniroot aggregate density hat qnorm model.matrix fivenum acf qt ecdf time approx embed frequency is.ts runif start ts optim quantile optimize dnorm dlnorm dexp dt t.test wilcox.test .preformat.ts var poly hclust as.dist
#' @importFrom utils globalVariables head tail combn flush.console
#' @importFrom xts to.monthly
#' @importFrom zoo as.yearmon index
#' @import data.table
#' @import doParallel
#' @import foreach
#' @rawNamespace import(Rcpp, except = LdFlags)
#' @import RcppParallel
#' @import rgl
#' @useDynLib NNS, .registration = TRUE



.onLoad <- function(libname = find.package("NNS"), pkgname = "NNS"){

  # CRAN Note avoidance

  utils::globalVariables(
    c("quadrant","quadrant.new","prior.quadrant",".","tmp.x","tmp.y","min_x_seg","max_x_seg","min_y_seg","max_y_seg",
      "mean_y_seg","mean_x_seg","sub.clpm",'sub.cupm','sub.dlpm','sub.dupm','weight','mean.x','mean.y',"upm","lpm","area",
      "Coefficient","X.Lower.Range","X.Upper.Range","y.hat","interval", "DISTANCES",
      "NNS.ID","max.x1","max.x2","min.x1","min.x2","counts",'old.counts',
      "Period","Coefficient.of.Variation","Variable.Coefficient.of.Variation", "Sum", "j","lpm","upm", "tau",
      "i.x","i.y","q_new","x.x","x.y","standard.errors",
      "detectCores","makeCluster", "makeForkCluster", "registerDoSEQ", "clusterExport", "frollmean", "shift",
      "%dopar%","foreach","stopCluster", "cl",
      "%do%", "k", "V1", "residuals", "nns_results", "bias_l", "bias_r",
      "bias", "conf.intervals", "conf.int.neg", "conf.int.pos", "pred.int", "lower.pred.int", "upper.pred.int",
      "estimates", "estimates.max", "estimates.min", "naive.first.grad", "naive.second.grad", "poly", "rise_1", "rise_2"
    ))

  requireNamespace("data.table")
  requireNamespace("doParallel")
  requireNamespace("foreach")
  requireNamespace("Rcpp")
  requireNamespace("RcppParallel")
  requireNamespace("rgl")
  

  .datatable.aware = TRUE
  
  options(datatable.verbose=FALSE)
  
  invisible(data.table::setDTthreads(0, throttle = NULL))
}

# Import calls and globalvariable calls

#' @importFrom grDevices adjustcolor rainbow rgb
#' @importFrom graphics abline boxplot legend lines par plot points segments text matplot title axis mtext barplot hist strwidth
#' @importFrom stats coef cor lm na.omit sd median complete.cases resid uniroot aggregate density hat qnorm model.matrix fivenum acf qt ecdf time approx embed frequency is.ts runif start ts optim quantile optimize dnorm dlnorm dexp dt .preformat.ts
#' @importFrom utils globalVariables head tail combn flush.console
#' @importFrom data.table data.table %chin% .I .N .SD := as.data.table fwrite is.data.table rbindlist set setcolorder setnames setorder as.IDate as.ITime setkey frollmean shift transpose
#' @importFrom dtw dtw
#' @importFrom meboot meboot
#' @importFrom tdigest tdigest tquantile
#' @importFrom Rfast colsums colmeans rowsums rowmeans
#' @importFrom caret upSample downSample
#' @importFrom plyr is.discrete
#' @importFrom zoo as.yearmon
#' @import doParallel
#' @import rgl
#' @import stringr
#' @import meboot
#' @import tdigest
#' @import data.table
#' @import dynlm
#' @import Quandl





.onLoad <- function(libname = find.package("NNS"), pkgname = "NNS"){

  # CRAN Note avoidance

  utils::globalVariables(
    c("quadrant","quadrant.new","prior.quadrant",".","tmp.x","tmp.y","min_x_seg","max_x_seg","min_y_seg","max_y_seg",
      "mean_y_seg","mean_x_seg","sub.clpm",'sub.cupm','sub.dlpm','sub.dupm','weight','mean.x','mean.y',"upm","lpm","area",
      "Coefficient","X.Lower.Range","X.Upper.Range","y.hat","interval", "DISTANCES",
      "NNS.ID","max.x1","max.x2","min.x1","min.x2","counts",'old.counts',
      "Period","Coefficient.of.Variation","Variable.Coefficient.of.Variation", "Sum", "j","lpm","upm", "tau",
      "i.x","i.y","q_new","x.x","x.y","standard.errors",
      "detectCores","makeCluster","registerDoSEQ","clusterExport", "frollmean", "shift",
      "%dopar%","foreach","stopCluster",
      "%do%", "k", "V1", "residuals", "nns_results", "bias_l", "bias_r",
      "tdigest", "tquantile"
    ))

  requireNamespace("data.table")
  requireNamespace("rgl")
  requireNamespace("doParallel")
  requireNamespace("stringr")
  requireNamespace("meboot")
  requireNamespace("tdigest")
  requireNamespace("Rfast")
  requireNamespace("dynlm")
  requireNamespace("Quandl")

  .datatable.aware = TRUE

  invisible()


}

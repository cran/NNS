#' NNS ANOVA
#'
#' Analysis of variance (ANOVA) based on lower partial moment CDFs for multiple variables.  Returns a degree of certainty the difference in sample means is zero, not a p-value.
#'
#' @param control a numeric vector, matrix or data frame.
#' @param treatment \code{NULL} (default) a numeric vector, matrix or data frame.
#' @param confidence.interval numeric [0, 1]; The confidence interval surrounding the \code{control} mean, defaults to \code{(confidence.interval = 0.95)}.
#' @param tails options: ("Left", "Right", "Both").  \code{tails = "Both"}(Default) Selects the tail of the distribution to determine effect size.
#' @param pairwise logical; \code{FALSE} (default) Returns pairwise certainty tests when set to \code{pairwise = TRUE}.
#' @param plot logical; \code{TRUE} (default) Returns the boxplot of all variables along with grand mean identification and confidence interval thereof.
#' @return Returns the following:
#' \itemize{
#' \item{\code{"Control Mean"}} \code{control} mean.
#' \item{\code{"Treatment Mean"}} \code{treatment} mean.
#' \item{\code{"Grand Mean"}} mean of means.
#' \item{\code{"Control CDF"}} CDF of the \code{control} from the grand mean.
#' \item{\code{"Treatment CDF"}} CDF of the \code{treatment} from the grand mean.
#' \item{\code{"Certainty"}} the certainty of the same population statistic.
#' \item{\code{"Lower Bound Effect"} and \code{"Upper Bound Effect"}} the effect size of the \code{treatment} for the specified confidence interval.
#' }
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{https://www.amazon.com/dp/1490523995/ref=cm_sw_su_dp}
#'
#' Viole, F. (2017) "Continuous CDFs and ANOVA with NNS"
#' \url{https://www.ssrn.com/abstract=3007373}
#'
#' @examples
#' ### Binary analysis and effect size
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.ANOVA(control = x, treatment = y)
#'
#' ### Two variable analysis with no control variable
#' A <- cbind(x, y)
#' NNS.ANOVA(A)
#'
#' ### Multiple variable analysis with no control variable
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100) ; z <- rnorm(100)
#' A <- cbind(x, y, z)
#' NNS.ANOVA(A)
#' @export


NNS.ANOVA <- function(control,
                     treatment,
                     confidence.interval = 0.95,
                     tails = "Both",
                     pairwise = FALSE,
                     plot = TRUE){

    tails <- tolower(tails)
    if(!any(tails%in%c("left","right","both"))){
        stop("Please select tails from 'left', 'right', or 'both'")
    }

    if(missing(treatment)){
        n <- ncol(control)
        if(is.null(n)){
            stop("supply both 'control' and 'treatment' or a matrix-like 'control'")
        }
        if(n == 1){
            stop("supply both 'control' and 'treatment' or a matrix-like 'control'")
        }
        if(n >= 2){
            if(any(class(control)=="tbl")) A <- as.data.frame(control)
        }
    } else {
        if(any(class(control)=="tbl")) control <- as.vector(unlist(control))
        if(any(class(treatment)=="tbl")) treatment <- as.vector(unlist(treatment))

        return(NNS.ANOVA.bin(control, treatment, confidence.interval = confidence.interval, plot = plot, tails = tails))
    }


    mean.of.means <- mean(colMeans(A))
    n <- ncol(A)
    if(!pairwise){
    #Continuous CDF for each variable from Mean of Means
        LPM_ratio <- sapply(1 : ncol(A), function(b) LPM.ratio(1, mean.of.means, A[ , b]))

        lower.25.target <- mean(sapply(1:n, function(i) LPM.VaR(.25,1,A[,i])))
        upper.25.target <- mean(sapply(1:n, function(i) UPM.VaR(.25,1,A[,i])))
        lower.125.target <- mean(sapply(1:n, function(i) LPM.VaR(.125,1,A[,i])))
        upper.125.target <- mean(sapply(1:n, function(i) UPM.VaR(.125,1,A[,i])))


        n <- ncol(A)
        raw.certainties <- list(n - 1)

        for(i in 1:(n - 1)){
          raw.certainties[[i]] <- sapply((i + 1) : n, function(b) NNS.ANOVA.bin(A[ , i],A[ , b],
                                                                                mean.of.means = mean.of.means,
                                                                                upper.25.target = upper.25.target,
                                                                                lower.25.target = lower.25.target,
                                                                                upper.125.target = upper.125.target,
                                                                                lower.125.target = lower.125.target,
                                                                                plot = FALSE)$Certainty)
        }


  #Certainty associated with samples
        NNS.ANOVA.rho <- mean(unlist(raw.certainties))

  #Graphs
        if(plot){
            boxplot(A, las = 2, xlab = "Means", ylab = "Variable", horizontal = TRUE, main = "NNS ANOVA", col = c('steelblue', rainbow(n - 1)))

  #For ANOVA Visualization
            abline(v = mean.of.means, col = "red", lwd = 4)
            mtext("Grand Mean", side = 3,col = "red", at = mean.of.means)
        }

        return(c("Certainty" = NNS.ANOVA.rho))

    } else {
          n <- ncol(A)
          raw.certainties <- list(n - 1)

          for(i in 1:(n - 1)){
              raw.certainties[[i]] <- sapply((i + 1) : n, function(b) NNS.ANOVA.bin(A[ , i],A[ , b], plot = FALSE)$Certainty)
          }

          certainties <- matrix(NA, n, n)
          certainties[lower.tri(certainties, diag = FALSE)] <- unlist(raw.certainties)
          diag(certainties) <- 1
          certainties <- pmax(certainties, t(certainties), na.rm = TRUE)

          colnames(certainties) <- colnames(A)
          rownames(certainties) <- colnames(A)

          if(plot){
              boxplot(A, las = 2, xlab = "Means", ylab = "Variable", horizontal = TRUE, main = "ANOVA", col = c('steelblue', rainbow(n - 1)))
              abline(v = mean.of.means, col = "red", lwd = 4)
              mtext("Grand Mean", side = 3,col = "red", at = mean.of.means)
          }

          return(certainties)
      }
}

## ----setup, include=FALSE, message=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(NNS)
library(data.table)
data.table::setDTthreads(2L)
options(mc.cores = 1)
Sys.setenv("OMP_THREAD_LIMIT" = 2)

## ----setup2, message=FALSE, warning = FALSE-----------------------------------
library(NNS)
library(data.table)
require(knitr)
require(rgl)

## -----------------------------------------------------------------------------
set.seed(123); x = rnorm(100)
ecdf(x)
P = ecdf(x)
P(0); P(1)

## ----message=FALSE------------------------------------------------------------
LPM.ratio(degree = 0, target = 0, variable = x); LPM.ratio(degree = 0, target = 1, variable = x)

## ----fig.align='center', fig.width=6, fig.height=6, echo = FALSE--------------
LPM.CDF = LPM.ratio(degree = 0, target = sort(x), variable = x)

plot(ecdf(x))
points(sort(x), LPM.CDF, col='red')
legend('left', legend = c('ecdf', 'LPM.ratio'), fill=c('black','red'), border=NA, bty='n')

## ----fig.align='center', fig.height=8, fig.width=8, echo=FALSE, warning=FALSE, message = FALSE, eval=FALSE----
#  zzz= rnorm(length(x), mean = 0, sd = 1)
#  norm_approx = pnorm(sort(zzz), mean=0, sd=1) #pnorm(sort(x),mean=-mean(x),sd=sd(x))
#  
#  plot(ecdf(x), main = "eCDF via LPM.ratio()", lwd = 4)
#  
#  
#  # Altering shape of distribution with LPM degree
#  for(i in c(0, 0.25, .5, 1, 2)){
#    idx <- which(i == c(0, 0.25, .5, 1, 2))
#    lines(sort(x), LPM.ratio(i, sort(x),x), col = rainbow(5, alpha = 1)[idx], lty = 1, lwd = 3)
#  }
#  
#   lines(sort(zzz), norm_approx ,col='black', lty = 3, lwd = 2)
#  
#  
#  legend("topleft",c("LPM.ratio(degree = 0)","LPM.ratio(degree = 0.25)","LPM.ratio(degree = 0.5)","LPM.ratio(degree = 1)","LPM.ratio(degree = 2)", "N(0,1) approximation"),
#         col = c(rainbow(5)[1:5], "black"), lwd = 3, lty = c(rep(1, 5), 3))

## ----fig.align='center', echo=FALSE, fig.width=10, fig.height=8, message=FALSE, warning=FALSE, eval=FALSE----
#  layout(matrix(c(1, 1, 1,1,1,
#                  2, 3, 4,5,6,
#                  2, 3, 4,5,6), nrow=5, byrow=FALSE),widths = c(2,rep(1,5)))
#  
#  
#  plot(ecdf(x), main = "eCDF via LPM.ratio()", lwd = 4)
#  
#  
#  # Altering shape of distribution with LPM degree
#  for(i in c(0, 0.25, .5, 1, 2)){
#    idx <- which(i == c(0, 0.25, .5, 1, 2))
#    lines(sort(x), LPM.ratio(i, sort(x),x), col = rainbow(5, alpha = 1)[idx], lty = 1, lwd = 3)
#  }
#  
#   lines(sort(zzz), norm_approx ,col='black', lty = 3, lwd = 2)
#  
#  
#  legend("topleft",c("LPM.ratio(degree = 0)","LPM.ratio(degree = 0.25)","LPM.ratio(degree = 0.5)","LPM.ratio(degree = 1)","LPM.ratio(degree = 2)", "N(0,1) approximation"),
#         col = c(rainbow(5)[1:5], "black"), lwd = 3, lty = c(rep(1, 5), 3))
#  
#  
#  
#  
#  y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
#  
#  plot(y$breaks,
#       c(y$counts,0), type = "s",
#      col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 0)", breaks = 15, xlab = "x", ylab = "freq")
#  hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), add = TRUE, col =  rainbow(5, alpha = .5)[1], breaks = 15)
#  
#  y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), border = NA, plot = FALSE, breaks = 15)
#  plot(y$breaks,
#       c(y$counts,0)
#       ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 0.25)", breaks = 15, xlab = "x", ylab = "freq")
#  hist(LPM.VaR(seq(0,1,length.out = 100), .25, x), border = rainbow(5)[2], add = TRUE, col =  rainbow(5, alpha = .5)[2], breaks = 15)
#  
#  y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
#  plot(y$breaks,
#       c(y$counts,0)
#       ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 0.5)", breaks = 15, xlab = "x", ylab = "freq")
#  hist(LPM.VaR(seq(0,1,length.out = 100), .5, x), border = rainbow(5)[3], add = TRUE, col =  rainbow(5, alpha = .5)[3], breaks = 15)
#  
#  y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
#  plot(y$breaks,
#       c(y$counts,0)
#       ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 1)", breaks = 15, xlab = "x", ylab = "freq")
#  hist(LPM.VaR(seq(0,1,length.out = 100), 1, x), border = rainbow(5)[4], add = TRUE, col =  rainbow(5, alpha = .5)[4], breaks = 15)
#  
#  y = hist(LPM.VaR(seq(0,1,length.out = 100), 0, x), plot = FALSE, breaks = 15)
#  plot(y$breaks,
#       c(y$counts,0)
#       ,type="s",col="black",lwd = 3, ylim = c(0,50), main = "Inverse CDF via LPM.VaR(degree 2)", breaks = 15, xlab = "x", ylab = "freq")
#  hist(LPM.VaR(seq(0,1,length.out = 100), 2, x), border = rainbow(5)[5], add = TRUE, col =  rainbow(5, alpha = .5)[5], breaks = 15)

## ----eval=FALSE---------------------------------------------------------------
#  degree.0.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 0, x = x)
#  degree.0.25.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 0.25, x = x)
#  degree.0.5.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 0.5, x = x)
#  degree.1.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 1, x = x)
#  degree.2.samples = LPM.VaR(percentile = seq(0, 1, length.out = 100), degree = 2, x = x)
#  
#  head(data.table::data.table(cbind("original x" = sort(x), degree.0.samples,
#                                                            degree.0.25.samples,
#                                                            degree.0.5.samples,
#                                                            degree.1.samples,
#                                                            degree.2.samples)), 10)
#  
#       original x degree.0.samples degree.0.25.samples degree.0.5.samples
#    1:  -2.309169        -2.309169           -2.309097         -2.3090915
#    2:  -1.966617        -1.966617           -1.941190         -1.6935509
#    3:  -1.686693        -1.686693           -1.599486         -1.4541494
#    4:  -1.548753        -1.548753           -1.382553         -1.2462731
#    5:  -1.265396        -1.265396           -1.250823         -1.1453748
#    6:  -1.265061        -1.265061           -1.176436         -1.0745440
#    7:  -1.220718        -1.220718           -1.119655         -1.0252742
#    8:  -1.138137        -1.138137           -1.067793         -0.9868693
#    9:  -1.123109        -1.123109           -1.026429         -0.9322105
#   10:  -1.071791        -1.071791           -1.014276         -0.8710942
#       degree.1.samples degree.2.samples
#    1:       -2.3091021       -2.3091170
#    2:       -1.4744653       -1.1614908
#    3:       -1.2159961       -0.9709972
#    4:       -1.0823023       -0.8610192
#    5:       -0.9968028       -0.7810300
#    6:       -0.9290505       -0.7169770
#    7:       -0.8666886       -0.6631888
#    8:       -0.8090433       -0.6170691
#    9:       -0.7556644       -0.5765608
#   10:       -0.7069835       -0.5403318

## ----fig.align='center', fig.width=8, fig.height=8, eval=FALSE----------------
#  boots = NNS.MC(x, reps = 1, lower_rho = -1, upper_rho = 1, by = .5)$replicates
#  reps = do.call(cbind, boots)
#  
#  plot(x, type = "l", lwd = 3, ylim = c(min(reps), max(reps)))
#  matplot(reps, type = "l", col = rainbow(length(boots)), add = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  sapply(boots, function(r) cor(r, x, method = "spearman"))
#  
#     rho = 1  rho = 0.5 rho = -0.5   rho = -1
#   1.0000000  0.4989619 -0.4984818 -0.9779778

## ----tgt_drift, fig.align='center', fig.width=8, fig.height=8, eval=FALSE-----
#  boots = NNS.MC(x, reps = 1, lower_rho = -1, upper_rho = 1, by = .5, target_drift = 0.05)$replicates
#  reps = do.call(cbind, boots)
#  
#  plot(x, type = "l", lwd = 3, ylim = c(min(reps), max(reps)))
#  matplot(reps, type = "l", col = rainbow(length(boots)), add = TRUE)

## ----multisim, eval=FALSE-----------------------------------------------------
#  set.seed(123)
#  x <- rnorm(1000); y <- rnorm(1000); z <- rnorm(1000)
#  
#  # Add variable x to original data to avoid total independence (example only)
#  original.data <- cbind(x, y, z, x)
#  
#  # Determine dependence structure
#  dep.structure <- apply(original.data, 2, function(x) LPM.ratio(degree = 1, target = x, variable = x))
#  
#  # Generate new data with different mean, sd and length (or distribution type)
#  new.data <- sapply(1:ncol(original.data), function(x) rnorm(nrow(original.data)*2, mean = 10, sd = 20))
#  
#  # Apply dependence structure to new data
#  new.dep.data <- sapply(1:ncol(original.data), function(x) LPM.VaR(percentile = dep.structure[,x], degree = 1, x = new.data[,x]))

## ----comparison, warning=FALSE, eval=FALSE------------------------------------
#  NNS.copula(original.data)
#  NNS.copula(new.dep.data)
#  
#  [1] 0.4353849
#  [1] 0.4357026

## ----eval=FALSE---------------------------------------------------------------
#  head(original.data)
#  head(new.dep.data)
#  
#                 x           y          z           x
#  [1,] -0.56047565 -0.99579872 -0.5116037 -0.56047565
#  [2,] -0.23017749 -1.03995504  0.2369379 -0.23017749
#  [3,]  1.55870831 -0.01798024 -0.5415892  1.55870831
#  [4,]  0.07050839 -0.13217513  1.2192276  0.07050839
#  [5,]  0.12928774 -2.54934277  0.1741359  0.12928774
#  [6,]  1.71506499  1.04057346 -0.6152683  1.71506499
#            [,1]       [,2]       [,3]      [,4]
#  [1,] -2.028109 -10.498044 -0.2090467 -1.682949
#  [2,]  4.608303 -11.390485 15.6213689  4.852534
#  [3,] 39.478741   8.836581 -0.8508203 40.585505
#  [4,] 10.683731   6.609255 36.0328589 10.877677
#  [5,] 11.866922 -47.955235 14.3111350 12.064633
#  [6,] 42.665726  29.639640 -2.4141874 43.797025

## ----eval=FALSE---------------------------------------------------------------
#  # Apply bootstrap to each variable
#  new.boot.dep.data = apply(original.data, 2, function(r) NNS.meboot(r, reps = 1, rho = .95))
#  
#  # Reformat into vectors
#  boot.ensemble.vectors = lapply(new.boot.dep.data, function(z) unlist(z["ensemble",]))
#  
#  # Create matrix from vectors
#  new.boot.dep.matrix = do.call(cbind, boot.ensemble.vectors)

## ----eval=FALSE---------------------------------------------------------------
#  for(i in 1:4) print(cor(new.boot.dep.matrix[,i], original.data[,i], method = "spearman"))
#  
#  [1] 0.9432899
#  [1] 0.9460947
#  [1] 0.9442031
#  [1] 0.9423242

## ----eval=FALSE---------------------------------------------------------------
#  NNS.copula(original.data)
#  NNS.copula(new.boot.dep.matrix)
#  
#  [1] 0.4353849
#  [1] 0.4263725

## ----eval=FALSE---------------------------------------------------------------
#  head(original.data)
#  head(new.boot.dep.matrix)
#  
#                 x           y          z           x
#  [1,] -0.56047565 -0.99579872 -0.5116037 -0.56047565
#  [2,] -0.23017749 -1.03995504  0.2369379 -0.23017749
#  [3,]  1.55870831 -0.01798024 -0.5415892  1.55870831
#  [4,]  0.07050839 -0.13217513  1.2192276  0.07050839
#  [5,]  0.12928774 -2.54934277  0.1741359  0.12928774
#  [6,]  1.71506499  1.04057346 -0.6152683  1.71506499
#                     x          y          z          x
#  ensemble1 -0.4268047 -0.7794553 -0.6364458 -0.4642642
#  ensemble2 -0.2965744 -1.0682197  0.3297265 -0.2531178
#  ensemble3  1.3302149  0.3054734 -0.4014515  1.4914884
#  ensemble4  0.2257378  0.3108846  1.0603892  0.1728540
#  ensemble5  0.4716743 -3.3344967 -0.1917697  0.4309379
#  ensemble6  1.3984978  1.1881374 -0.5295386  1.5326055

## ----threads, echo = FALSE----------------------------------------------------
Sys.setenv("OMP_THREAD_LIMIT" = "")


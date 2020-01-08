## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----setup2,message=FALSE-----------------------------------------------------
require(NNS)
require(knitr)
require(rgl)
require(data.table)
require(dtw)

## ----linear,fig.width=5,fig.height=3,fig.align = "center"---------------------
x = seq(0, 3, .01) ; y = 2 * x

cor(x, y)
NNS.dep(x, y, print.map = TRUE, order = 3)

## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center"------------------
x=seq(0, 3, .01) ; y = x ^ 10

cor(x, y)
NNS.dep(x, y, print.map = TRUE)

## ----dependence,fig.width=5,fig.height=3,fig.align = "center"-----------------
set.seed(123)
df <- data.frame(x = runif(10000, -1, 1), y = runif(10000, -1, 1))
df <- subset(df, (x ^ 2 + y ^ 2 <= 1 & x ^ 2 + y ^ 2 >= 0.95))
NNS.dep(df$x, df$y, print.map = TRUE)

## ----multi--------------------------------------------------------------------
set.seed(123)
x <- rnorm(1000); y <- rnorm(1000); z <- rnorm(1000)
NNS.dep.hd(cbind(x, y, z), plot = TRUE, independence.overlay = TRUE)

## ----permutations-------------------------------------------------------------
## p-values for [NNS.dep]
x <- seq(-5, 5, .1); y <- x^2 + rnorm(length(x))

nns_cor_dep <- NNS.dep(x, y, print.map = TRUE)
nns_cor_dep

## Create permutations of y
y_p <- replicate(100, sample.int(length(y)))

## Generate new correlation and dependence measures on each new permutation of y
nns.mc <- apply(y_p, 2, function(g) NNS.dep(x, y[g]))

## Store results
cors <- unlist(lapply(nns.mc, "[[", 1))
deps <- unlist(lapply(nns.mc, "[[", 2))

## View results
hist(cors)
abline(v = LPM.VaR(.975,0, cors), col = 'red')
abline(v = UPM.VaR(.975,0, cors), col = 'red')


## Left tailed correlation p-value
cor_p_value <- LPM(0, nns_cor_dep$Correlation, cors)
cor_p_value

## Right tailed correlation p-value
cor_p_value <- UPM(0, nns_cor_dep$Correlation, cors)
cor_p_value

## Confidence Intervals
## For 95th percentile VaR (both-tails) see [LPM.VaR] and [UPM.VaR]
## Lower CI
LPM.VaR(.975, 0, cors)
## Upper CI
UPM.VaR(.975, 0, cors)


hist(deps)
abline(v = LPM.VaR(.975,0, deps), col = 'red')
abline(v = UPM.VaR(.975,0, deps), col = 'red')


## Left tailed dependence p-value
dep_p_value <- LPM(0, nns_cor_dep$Dependence, deps)
dep_p_value

## Right tailed dependence p-value
dep_p_value <- UPM(0, nns_cor_dep$Dependence, deps)
dep_p_value

## Confidence Intervals
## For 95th percentile VaR (both-tails) see [LPM.VaR] and [UPM.VaR]
## Lower CI
LPM.VaR(.975, 0, deps)
## Upper CI
UPM.VaR(.975, 0, deps)


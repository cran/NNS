## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----setup2, message=FALSE-----------------------------------------------
require(NNS)
require(knitr)
require(rgl)
require(data.table)

## ----rhs-----------------------------------------------------------------
NNS.reg(iris[,1:4], iris[,5], residual.plot = FALSE)$rhs.partitions

## ----NNSBOOST,fig.align = "center", fig.height = 8,fig.width=6.5---------
set.seed(123)
test.set = sample(150,10)

a = NNS.boost(iris[-test.set, 1:4], iris[-test.set, 5],
              IVs.test = iris[test.set, 1:4],
              epochs = 100, learner.trials = 100, status = FALSE)

mean(round(a)==as.numeric(iris[test.set,5]))


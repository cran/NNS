## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----setup2, message=FALSE-----------------------------------------------
require(NNS)
require(knitr)
require(rgl)
require(data.table)

## ----linear,fig.width=5,fig.height=3,fig.align = "center", warning=FALSE----
nns = NNS.ARMA(AirPassengers, h = 44, training.set = 100, method = "lin", plot = TRUE, seasonal.factor = 12, seasonal.plot = FALSE, ncores = 1)

sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----nonlinear,fig.width=5,fig.height=3,fig.align = "center"-------------
nns = NNS.ARMA(AirPassengers, h = 44, training.set = 100, method = "nonlin", plot = FALSE, seasonal.factor = 12, seasonal.plot = FALSE, ncores = 1)

sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----seasonal test,fig.width=5,fig.height=3,fig.align = "center"---------
seas = t(sapply(1 : 25, function(i) c(i, sqrt( mean( (NNS.ARMA(AirPassengers, h = 44, training.set = 100, method = "lin", seasonal.factor = i, plot=FALSE, ncores = 1) - tail(AirPassengers, 44)) ^ 2) ) ) ) )

colnames(seas) = c("Period", "RMSE")
seas

## ----best fit------------------------------------------------------------
a = seas[which.min(seas[ , 2]), 1]

## ----best nonlinear,fig.width=5,fig.height=3,fig.align = "center"--------
nns = NNS.ARMA(AirPassengers, h = 44, training.set = 100, method = "nonlin", seasonal.factor = a, plot = TRUE, seasonal.plot = FALSE, ncores = 1)

sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----modulo--------------------------------------------------------------
NNS.seas(AirPassengers, modulo = 12, plot = FALSE)

## ----best optim,fig.width=5,fig.height=3,fig.align = "center"------------
nns.optimal = NNS.ARMA.optim(AirPassengers[1:100], 
                             training.set = 88, 
                             seasonal.factor = seq(12, 24, 6),
                             obj.fn = expression( sqrt(mean((predicted - actual)^2)) ), ncores = 1, subcores = 1)

nns.optimal

## ----best optim2,fig.width=5,fig.height=3,fig.align = "center"-----------
nns = NNS.ARMA(AirPassengers, training.set = 100, h = 44, seasonal.factor = nns.optimal$periods, weights = nns.optimal$weights, method = nns.optimal$method, plot = TRUE, seasonal.plot = FALSE, ncores = 1)

sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----best optim3,fig.width=5,fig.height=3,fig.align = "center"-----------
sqrt(mean((nns+nns.optimal$bias.shift - tail(AirPassengers, 44)) ^ 2))

## ----neg-----------------------------------------------------------------
nns <- pmax(0, nns+nns.optimal$bias.shift)
sqrt(mean((nns - tail(AirPassengers, 44)) ^ 2))

## ----extension,results='hide',fig.width=5,fig.height=3,fig.align = "center"----
NNS.ARMA(AirPassengers, h = 50, seasonal.factor = nns.optimal$periods, method  = nns.optimal$method, weights = nns.optimal$weights, plot = TRUE, seasonal.plot = FALSE, ncores = 1) + nns.optimal$bias.shift


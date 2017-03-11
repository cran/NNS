## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
require(NNS)
require(knitr)
require(rgl)

## ----linear--------------------------------------------------------------
x=seq(0,3,.01); y=2*x

cor(x,y)
NNS.dep(x,y,print.map = T)

## ----nonlinear-----------------------------------------------------------
x=seq(0,3,.01); y=x^10

cor(x,y)
NNS.dep(x,y,print.map = T)

## ----dependence----------------------------------------------------------
set.seed(123)
df<- data.frame(x=runif(10000,-1,1),y=runif(10000,-1,1))
df<- subset(df, (x^2 + y^2 <= 1 & x^2 + y^2 >= 0.95))
NNS.dep(df$x,df$y,print.map = T)

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Stack.R
\name{NNS.stack}
\alias{NNS.stack}
\title{NNS Stack}
\usage{
NNS.stack(
  IVs.train,
  DV.train,
  IVs.test = NULL,
  type = NULL,
  obj.fn = expression(sum((predicted - actual)^2)),
  objective = "min",
  optimize.threshold = TRUE,
  dist = "L2",
  CV.size = NULL,
  balance = FALSE,
  ts.test = NULL,
  folds = 5,
  order = NULL,
  norm = NULL,
  method = c(1, 2),
  stack = TRUE,
  dim.red.method = "cor",
  pred.int = NULL,
  status = TRUE,
  ncores = NULL
)
}
\arguments{
\item{IVs.train}{a vector, matrix or data frame of variables of numeric or factor data types.}

\item{DV.train}{a numeric or factor vector with compatible dimensions to \code{(IVs.train)}.}

\item{IVs.test}{a vector, matrix or data frame of variables of numeric or factor data types with compatible dimensions to \code{(IVs.train)}.  If NULL, will use \code{(IVs.train)} as default.}

\item{type}{\code{NULL} (default).  To perform a classification of discrete integer classes from factor target variable \code{(DV.train)} with a base category of 1, set to \code{(type = "CLASS")}, else for continuous \code{(DV.train)} set to \code{(type = NULL)}.   Like a logistic regression, this setting is not necessary for target variable of two classes e.g. [0, 1].}

\item{obj.fn}{expression; \code{expression(sum((predicted - actual)^2))} (default) Sum of squared errors is the default objective function.  Any \code{expression()} using the specific terms \code{predicted} and \code{actual} can be used.}

\item{objective}{options: ("min", "max") \code{"min"} (default) Select whether to minimize or maximize the objective function \code{obj.fn}.}

\item{optimize.threshold}{logical; \code{TRUE} (default) Will optimize the probability threshold value for rounding in classification problems.  If \code{FALSE}, returns 0.5.}

\item{dist}{options:("L1", "L2", "DTW", "FACTOR") the method of distance calculation; Selects the distance calculation used. \code{dist = "L2"} (default) selects the Euclidean distance and \code{(dist = "L1")} selects the Manhattan distance; \code{(dist = "DTW")} selects the dynamic time warping distance; \code{(dist = "FACTOR")} uses a frequency.}

\item{CV.size}{numeric [0, 1]; \code{NULL} (default) Sets the cross-validation size if \code{(IVs.test = NULL)}.  Defaults to a random value between 0.2 and 0.33 for a random sampling of the training set.}

\item{balance}{logical; \code{FALSE} (default) Uses both up and down sampling to balance the classes.  \code{type="CLASS"} required.}

\item{ts.test}{integer; NULL (default) Sets the length of the test set for time-series data; typically \code{2*h} parameter value from \link{NNS.ARMA} or double known periods to forecast.}

\item{folds}{integer; \code{folds = 5} (default) Select the number of cross-validation folds.}

\item{order}{options: (integer, "max", NULL); \code{NULL} (default) Sets the order for \link{NNS.reg}, where \code{(order = "max")} is the k-nearest neighbors equivalent, which is suggested for mixed continuous and discrete (unordered, ordered) data.}

\item{norm}{options: ("std", "NNS", NULL); \code{NULL} (default) 3 settings offered: \code{NULL}, \code{"std"}, and \code{"NNS"}.  Selects the \code{norm} parameter in \link{NNS.reg}.}

\item{method}{numeric options: (1, 2); Select the NNS method to include in stack.  \code{(method = 1)} selects \link{NNS.reg}; \code{(method = 2)} selects \link{NNS.reg} dimension reduction regression.  Defaults to \code{method = c(1, 2)}, which will reduce the dimension first, then find the optimal \code{n.best}.}

\item{stack}{logical; \code{TRUE} (default) Uses dimension reduction output in \code{n.best} optimization, otherwise performs both analyses independently.}

\item{dim.red.method}{options: ("cor", "NNS.dep", "NNS.caus", "equal", "all") method for determining synthetic X* coefficients.  \code{(dim.red.method = "cor")} uses standard linear correlation for weights.  \code{(dim.red.method = "NNS.dep")} (default) uses \link{NNS.dep} for nonlinear dependence weights, while \code{(dim.red.method = "NNS.caus")} uses \link{NNS.caus} for causal weights.  \code{(dim.red.method = "all")} averages all methods for further feature engineering.}

\item{pred.int}{numeric [0,1]; \code{NULL} (default) Returns the associated prediction intervals with each \code{method}.}

\item{status}{logical; \code{TRUE} (default) Prints status update message in console.}

\item{ncores}{integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.reg}. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.}
}
\value{
Returns a vector of fitted values for the dependent variable test set for all models.
\itemize{
\item{\code{"NNS.reg.n.best"}} returns the optimum \code{"n.best"} parameter for the \link{NNS.reg} multivariate regression.  \code{"SSE.reg"} returns the SSE for the \link{NNS.reg} multivariate regression.
\item{\code{"OBJfn.reg"}} returns the \code{obj.fn} for the \link{NNS.reg} regression.
\item{\code{"NNS.dim.red.threshold"}} returns the optimum \code{"threshold"} from the \link{NNS.reg} dimension reduction regression.
\item{\code{"OBJfn.dim.red"}} returns the \code{obj.fn} for the \link{NNS.reg} dimension reduction regression.
\item{\code{"probability.threshold"}} returns the optimum probability threshold for classification, else 0.5 when set to \code{FALSE}.
\item{\code{"reg"}} returns \link{NNS.reg} output.
\item{\code{"reg.pred.int"}} returns the prediction intervals for the regression output.
\item{\code{"dim.red"}} returns \link{NNS.reg} dimension reduction regression output.
\item{\code{"dim.red.pred.int"}} returns the prediction intervals for the dimension reduction regression output.
\item{\code{"stack"}} returns the output of the stacked model.
\item{\code{"pred.int"}} returns the prediction intervals for the stacked model.
}
}
\description{
Prediction model using the predictions of the NNS base models \link{NNS.reg} as features (i.e. meta-features) for the stacked model.
}
\note{
\itemize{
\item Incorporate any objective function from external packages (such as \code{Metrics::mape}) via \code{NNS.stack(..., obj.fn = expression(Metrics::mape(actual, predicted)), objective = "min")}

\item Like a logistic regression, the \code{(type = "CLASS")} setting is not necessary for target variable of two classes e.g. [0, 1].  The response variable base category should be 1 for multiple class problems.

\item Missing data should be handled prior as well using \link{na.omit} or \link{complete.cases} on the full dataset.
}

If error received:

\code{"Error in is.data.frame(x) : object 'RP' not found"}

reduce the \code{CV.size}.
}
\examples{
 ## Using 'iris' dataset where test set [IVs.test] is 'iris' rows 141:150.
 \dontrun{
 NNS.stack(iris[1:140, 1:4], iris[1:140, 5], IVs.test = iris[141:150, 1:4], type = "CLASS")

 ## Using 'iris' dataset to determine [n.best] and [threshold] with no test set.
 NNS.stack(iris[ , 1:4], iris[ , 5], type = "CLASS")

 ## Selecting NNS.reg and dimension reduction techniques.
 NNS.stack(iris[1:140, 1:4], iris[1:140, 5], iris[141:150, 1:4], method = c(1, 2), type = "CLASS")
 }
}
\references{
Viole, F. (2016) "Classification Using NNS Clustering Analysis"  \doi{10.2139/ssrn.2864711}
}
\author{
Fred Viole, OVVO Financial Systems
}

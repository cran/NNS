% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Normalization.R
\name{NNS.norm}
\alias{NNS.norm}
\title{NNS Normalization}
\usage{
NNS.norm(X, linear = FALSE, chart.type = NULL, location = "topleft")
}
\arguments{
\item{X}{a numeric matrix or data frame, or a list.}

\item{linear}{logical; \code{FALSE} (default) Performs a linear scaling normalization, resulting in equal means for all variables.}

\item{chart.type}{options: ("l", "b"); \code{NULL} (default).  Set \code{(chart.type = "l")} for line,
\code{(chart.type = "b")} for boxplot.}

\item{location}{Sets the legend location within the plot, per the \code{x} and \code{y} co-ordinates used in base graphics \link{legend}.}
}
\value{
Returns a \link{data.frame} of normalized values.
}
\description{
Normalizes a matrix of variables based on nonlinear scaling normalization method.
}
\note{
Unequal vectors provided in a list will only generate \code{linear=TRUE} normalized values.
}
\examples{
\dontrun{
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100)
A <- cbind(x, y)
NNS.norm(A)

### Normalize list of unequal vector lengths

vec1 <- c(1, 2, 3, 4, 5, 6, 7)
vec2 <- c(10, 20, 30, 40, 50, 60)
vec3 <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)

vec_list <- list(vec1, vec2, vec3)
NNS.norm(vec_list)
}
}
\references{
Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
}
\author{
Fred Viole, OVVO Financial Systems
}

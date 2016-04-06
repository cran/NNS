#' VN Normalization
#'
#' Normalizes a matrix of variables based on nonlinear scaling normalization method.
#' @param A Matrix of variables.
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' A<-cbind(x,y)
#' \dontrun{VN.norm(A)}
#' @export

VN.norm <- function(A) {
  m  <- colMeans(A)
  RG <- m %o% (1/m)
  scales <- colMeans(RG * abs(cor(A)))
  A_Normalized <- t(t(A) * scales)

  n <- ncol(A)
  i <- seq_len(n)
  labels <- c(sprintf("Variable%i", i),
              sprintf("Variable%i_Normalized", i))
  boxplot(cbind(A, A_Normalized),
          las = 2, names = labels,
          col = c(rep("grey", n), rainbow(n)))


  return(A_Normalized)

}

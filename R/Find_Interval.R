#' Find Interval
#'
#' Internal Function Only.
#' Removes the sort check from the base function \link{findInterval}
#' @param x Variable / Vector
#' @param v Vector
#' @return Returns interval number




findInterval2 <- function(x,v) {
  n = length(v)
  if (x<=v[1])
    return (1)
  if (x>=v[n])
    return (n+1)
  i=1
  k=n
  while({j = (k-i) %/% 2 + i; !(v[j] <= x && x < v[j+1])}) {
    if (x < v[j])
      k = j
    else
      i = j+1
  }
  return (j)
}



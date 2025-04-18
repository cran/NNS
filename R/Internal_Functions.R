### Continuous Mode of a distribution
mode <- function(x) NNS.mode(x, discrete = FALSE, multi = FALSE)


### Classification Mode of a distribution
mode_class <- function(x) NNS.mode(x, discrete = TRUE, multi = FALSE)


### Gravity of a distribution
gravity <- function(x) NNS.gravity(x, discrete = FALSE)

gravity_class <- function(x) NNS.gravity(x, discrete = TRUE)


### Factor to dummy variable
factor_2_dummy <- function(x){
  x <- unlist(x)
  if(is.factor(x) && length(unique(x)) > 1){
    output <- model.matrix(~(x) -1, x)[,-1]
  } else {
    output <- as.numeric(x)
  }
  output
}

### Factor to dummy variable FULL RANK
factor_2_dummy_FR <- function(x){
  x <- unlist(x)
  if(is.factor(x) && length(unique(x)) > 1){
    output <- model.matrix(~(x) -1, x)
  } else {
    output <- as.numeric(x)
  }
  output
}

### Generator for 1:length(lag) vectors in NNS.ARMA
generate.vectors <- function(x, l) {
  Component.series <- lapply(l, function(lag) {
    rev.series <- rev(x[seq(length(x) + 1, 1, -lag)])
    rev.series[!is.na(rev.series)]
  })
  
  Component.index <- lapply(Component.series, seq_along)
  
  return(list(Component.index = Component.index, Component.series = Component.series))
}


generate.lin.vectors <- function(x, l, h = 1) {
  original.index <- seq_along(x)
  augmented.index <- c(original.index, tail(original.index,1) + (1:h))
  max_fcast <- min(h, l)  
  # Generate lagged components by applying lag across the augmented index
  Component.series <- lapply(1:max_fcast, function(i) {
    start.index <- length(x) + i  # Shift by 1 each time
    rev.series <- rev(x[seq(start.index, 1, -l)])  # Reverse the sequence with a lag of 12
    rev.series[!is.na(rev.series)]  # Remove any NA values
  })
  
  Component.index <- lapply(Component.series, seq_along)
  
  # Initialize forecast.index and forecast.values
  forecast.index <- vector("list", length(l))
  forecast.values <- vector("list", length(l))
  
  # Generate forecast.index for 1:h
  max_fcast <- min(h, l)
  forecast.index <- create_recycled_list(1:h, max_fcast)
  forecast.index <- forecast.index[!sapply(forecast.index, is.null)]
  
  
  # Generate forecast values based on the last value in Component.index
  forecast.values.raw <- lapply(1:h, function(i) {
    # Recycle the index if h > l
    recycled_index <- (i - 1) %% l + 1  # Cycle through 1 to l
    
    # Get the last value from the corresponding Component.index
    last_value <- tail(Component.index[[recycled_index]], 1)
    
    # Calculate the forecast increment
    forecast_increment <- ceiling(i / l)
    
    # Generate the forecast value
    forecast_value <- last_value + forecast_increment
    return(forecast_value)
  })

  forecast.values <- create_recycled_list(unlist(forecast.values.raw), l)
  forecast.values <- forecast.values[!sapply(forecast.values, is.null)]
  
  return(list(
    Component.index = Component.index,
    Component.series = Component.series,
    forecast.values = forecast.values,
    forecast.index = forecast.index
  ))
}

create_recycled_list <- function(values, list_length) {
  # Initialize an empty list to store the values
  result <- vector("list", list_length)
  
  # Recycle the values by repeating them across the list length
  for (i in seq_along(values)) {
    index <- (i - 1) %% list_length + 1
    result[[index]] <- c(result[[index]], values[i])
  }
  
  return(result)
}


### Weight and lag function for seasonality in NNS.ARMA
ARMA.seas.weighting <- function(sf,mat){
  M <- mat
  n <- ncol(M)
  if(is.null(n)){
    return(list(lag = M[1], Weights = 1))
  }

  if(n == 1){
    return(list(lag = 1, Weights = 1))
  }

  if(n > 1){
    if(sf){
      lag <- M$all.periods$Period[1]
      Weights <- 1
      return(list(lag = lag, Weights = Weights))
    }

    # Determine lag from seasonality test
    if(!sf){
      lag <- na.omit(unlist(M$Period))
      Observation.weighting <- (1 / sqrt(lag))
      if(is.na(M$Coefficient.of.Variation)  && length(M$Coefficient.of.Variation)==1){
        Lag.weighting <- 1
      } else {
        Lag.weighting <- (unlist(M$Variable.Coefficient.of.Variation) - unlist(M$Coefficient.of.Variation))
      }
      Weights <- (Lag.weighting * Observation.weighting) / sum(Lag.weighting * Observation.weighting)
      return(list(lag = lag, Weights = Weights))
    }
  }
}


### Lag matrix generator for NNS.VAR
### Vector of tau for single different tau per variables tau = c(1, 4)
### List of tau vectors for multiple different tau per variables tau = list(c(1,2,3), c(4,5,6))
lag.mtx <- function(x, tau){
  colheads <- NULL
  
  max_tau <- max(unlist(tau))
  
  if(is.null(dim(x)[2])) {
    colheads <- noquote(as.character(deparse(substitute(x))))
    x <- t(t(x))
  }
  
  j.vectors <- vector(mode = "list", ncol(x))
  
  for(j in 1:ncol(x)){
    if(is.null(colheads)){
      colheads <- colnames(x)[j]
      
      colheads <- noquote(as.character(deparse(substitute(colheads))))
    }
    
    x.vectors <- vector(mode = "list")
    heads <- paste0(colheads, "_tau_")
    heads <- gsub('"', '' ,heads)
    
    for (i in 0:max_tau){
      x.vectors[[paste(heads, i, sep = "")]] <- numeric(0L)
      start <- max_tau - i + 1
      end <- length(x[,j]) - i
      x.vectors[[i + 1]] <- x[start : end, j]
    }
    
    j.vectors[[j]] <- do.call(cbind, x.vectors)
    colheads <- NULL
  }
  mtx <- as.data.frame(do.call(cbind, j.vectors))
  

  if(length(unlist(tau)) > 1){
    relevant_lags <- lapply(1:length(tau), function(i) c((i-1)*max_tau + i, (i-1)*max_tau + unlist(tau[[i]]) + i))

    relevant_lags <- sort(unlist(relevant_lags))
    mtx <- mtx[ , relevant_lags]
  }
  
  vars <- which(grepl("tau_0", colnames(mtx)))
  
  everything_else <- seq_len(dim(mtx)[2])[-vars]
  mtx <- mtx[,c(vars, everything_else)]
  
  return(mtx)
}




### Refactored meboot::meboot.part function
NNS.meboot.part <- function(x, n, z, xmin, xmax, desintxb, reachbnd)
{
  # Generate random numbers from the [0,1] uniform interval
  p <- runif(n, min=0, max=1)

  q <- quantile(x, p)

  ref1 <- which(p <= (1/n))
  if(length(ref1) > 0){
    qq <- approx(c(0,1/n), c(xmin,z[1]), p[ref1])$y
    q[ref1] <- qq
    if(!reachbnd)  q[ref1] <- qq + desintxb[1]-0.5*(z[1]+xmin)
  }

  ref4 <- which(p == ((n-1)/n))
  if(length(ref4) > 0)
    q[ref4] <- z[n-1]

  ref5 <- which(p > ((n-1)/n))
  if(length(ref5) > 0){
    # Right tail proportion p[i]
    qq <- approx(c((n-1)/n,1), c(z[n-1],xmax), p[ref5])$y
    q[ref5] <- qq   # this implicitly shifts xmax for algorithm
    if(!reachbnd)  q[ref5] <- qq + desintxb[n]-0.5*(z[n-1]+xmax)
    # such that the algorithm gives xmax when p[i]=1
    # this is the meaning of reaching the bounds xmax and xmin
  }

  q

}

### Refactored meboot::expand.sd function
NNS.meboot.expand.sd <- function(x, ensemble, fiv=5){
  sdx <- if (is.null(ncol(x))) sd(x) else apply(x, 2, sd)

  sdf <- c(sdx, apply(ensemble, 2, sd))

  sdfa <- sdf/sdf[1]  # ratio of actual sd to that of original data
  sdfd <- sdf[1]/sdf  # ratio of desired sd to actual sd

  # expansion is needed since some of these are <1 due to attenuation
  mx <- 1+(fiv/100)
  # following are expansion factors
  id <- which(sdfa < 1)
  if (length(id) > 0) sdfa[id] <- runif(n=length(id), min=1, max=mx)

  sdfdXsdfa <- sdfd[-1]*sdfa[-1]

  id <- which(floor(sdfdXsdfa) > 0)

  if (length(id) > 0) {
    if(length(id) > 1) ensemble[,id] <- ensemble[,id] %*% diag(sdfdXsdfa[id]) else ensemble[,id] <- ensemble[,id] * sdfdXsdfa[id]
  }

  if(is.ts(x)) ensemble <- ts(ensemble, frequency=frequency(x), start=start(x))


  ensemble
}


# Refactored force.clt function from meboot
force.clt <- function(x, ensemble)
{
  n <- nrow(ensemble)
  bigj <- ncol(ensemble)
  
  gm <- mean(x)  # desired grand mean
  
  s <- if (is.null(ncol(x)))
    sd(x) else apply(x, 2, sd)
  smean <- s/sqrt(bigj)  # desired standard deviation of means by CLT
  xbar <- apply(ensemble, 2, mean)
  sortxbar <- sort(xbar)  
  oo <- order(xbar)
  
  # now spread out the means as if they were from normal density
  # smallest mean should equal the normal quantile at prob=1/bigj+1
  # smallest second mean= the normal quantile at prob=2/bigj+1
  # . . .
  # LAST mean= the normal quantile at prob=bigj/(bigj+1)
  
  newbar <- gm + qnorm(1:bigj/(bigj+1))*smean
  
  # the above adjustement of means will change their sd
  # CLT says that their sd should be s / sqrt(n)
  # so we have recenter and rescale these revised means called newbar
  
  scn <- scale(newbar)  # first scale them to have zero mean and unit sd
  newm <- scn*smean+gm  # this forces the mean to be gm and sd=s / sqrt(n)
  
  meanfix <- as.numeric(newm - sortxbar)
  
  # we have lost the original order in sorting, need to go back 
  out <- ensemble
  for(i in 1:bigj)
    out[,oo[i]] <- ensemble[,oo[i]] + meanfix[i]
  
  if(any(is.ts(ensemble))){
    out <- ts(out, frequency=frequency(ensemble), start=start(ensemble))
    dimnames(out)[[2]] <- dimnames(ensemble)[[2]]
  }
  out
}


is.fcl <- function(x) is.factor(x) || is.character(x) || is.logical(x)

is.discrete <- function(x) sum(as.numeric(x)%%1)==0


### upSample / downSample to avoid dependencies
downSample <- function(x, y, list = FALSE, yname = "Class") {
  # Ensure x is a data frame
  if (!is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = TRUE)
  }
  # Ensure y is a factor
  if (!is.factor(y)) {
    warning(
      "Down-sampling requires a factor variable as the response. The original data was returned."
    )
    return(list(x = x, y = y))
  }
  
  # Determine the minimum class size
  minClass <- min(table(y))
  
  # Create an empty list to store sampled data
  sampled_data <- vector("list", length(unique(y)))
  names(sampled_data) <- unique(y)
  
  # Down-sample each class
  for (class in names(sampled_data)) {
    class_indices <- which(y == class)
    sampled_indices <- sample(class_indices, minClass)
    sampled_data[[class]] <- x[sampled_indices, , drop = FALSE]
  }
  
  # Combine the down-sampled data
  x <- do.call(rbind, sampled_data)
  
  # Extract the outcome and remove it from x
  y <- factor(rep(names(sampled_data), each = minClass))
  
  # Prepare the output
  if (list) {
    out <- list(x = x, y = y)
  } else {
    out <- cbind(x, y)
    colnames(out)[ncol(out)] <- yname
  }
  
  return(out)
}



upSample <- function(x, y, list = FALSE, yname = "Class") {
  # Ensure x is a data frame
  if (!is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = TRUE)
  }
  # Ensure y is a factor
  if (!is.factor(y)) {
    warning(
      "Up-sampling requires a factor variable as the response. The original data was returned."
    )
    return(list(x = x, y = y))
  }
  
  # Determine the maximum class size
  maxClass <- max(table(y))
  
  # Create an empty list to store sampled data
  sampled_data <- vector("list", length(unique(y)))
  names(sampled_data) <- unique(y)
  
  # Up-sample each class
  for (class in names(sampled_data)) {
    class_indices <- which(y == class)
    class_data <- x[class_indices, , drop = FALSE]
    if (nrow(class_data) < maxClass) {
      extra_indices <- sample(seq_len(nrow(class_data)), size = maxClass - nrow(class_data), replace = TRUE)
      sampled_data[[class]] <- rbind(class_data, class_data[extra_indices, , drop = FALSE])
    } else {
      sampled_data[[class]] <- class_data
    }
  }
  
  # Combine the up-sampled data
  x <- do.call(rbind, sampled_data)
  
  # Extract the outcome and remove it from x
  y <- factor(rep(names(sampled_data), each = maxClass))
  
  # Prepare the output
  if (list) {
    out <- list(x = x, y = y)
  } else {
    out <- cbind(x, y)
    colnames(out)[ncol(out)] <- yname
  }
  
  return(out)
}

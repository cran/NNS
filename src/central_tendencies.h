// central_tendencies.h
#ifndef CENTRAL_TENDENCIES_H
#define CENTRAL_TENDENCIES_H

#include <Rcpp.h>
using namespace Rcpp;

// Declare the functions without default values
double NNS_gravity_cpp(NumericVector x, bool discrete);
NumericVector NNS_mode_cpp(NumericVector x, bool discrete, bool multi);

#endif


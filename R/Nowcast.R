#' NNS Nowcast
#'
#' Wrapper function for NNS nowcasting method using the nonparametric vector autoregression \link{NNS.VAR}, and Federal Reserve Nowcasting variables.
#'
#' @param h integer; \code{(h = 1)} (default) Number of periods to forecast. \code{(h = 0)} will return just the interpolated and extrapolated values up to the current month.
#' @param additional.regressors character; \code{NULL} (default) add more regressors to the base model.  The format must utilize the \code{\link[quantmod]{getSymbols}} format for FRED data, else specify the source.
#' @param additional.sources character; \code{NULL} (default) specify the \code{source} argument per \code{\link[quantmod]{getSymbols}} for each \code{additional.regressors} specified.
#' @param naive.weights logical; \code{TRUE} Equal weights applied to univariate and multivariate outputs in ensemble.  \code{FALSE} (default) will apply weights based on the number of relevant variables detected. 
#' @param specific.regressors integer; \code{NULL} (default) Select individual regressors from the base model per Viole (2020) listed in the \code{Note} below.
#' @param start.date character; \code{"2000-01-03"} (default) Starting date for all data series download.
#' @param keep.data logical; \code{FALSE} (default) Keeps downloaded variables in a new environment \code{NNSdata}.
#' @param status logical; \code{TRUE} (default) Prints status update message in console.
#' @param ncores integer; value specifying the number of cores to be used in the parallelized subroutine \link{NNS.ARMA.optim}. If NULL (default), the number of cores to be used is equal to the number of cores of the machine - 1.
#'
#' @note Specific regressors include:
#' \enumerate{
#'   \item \code{PAYEMS} -- Payroll Employment
#'   \item \code{JTSJOL} -- Job Openings
#'   \item \code{CPIAUCSL} -- Consumer Price Index
#'   \item \code{DGORDER} -- Durable Goods Orders
#'   \item \code{RSAFS} -- Retail Sales
#'   \item \code{UNRATE} -- Unemployment Rate
#'   \item \code{HOUST} -- Housing Starts
#'   \item \code{INDPRO} -- Industrial Production
#'   \item \code{DSPIC96} -- Personal Income
#'   \item \code{BOPTEXP} -- Exports
#'   \item \code{BOPTIMP} -- Imports
#'   \item \code{TTLCONS} -- Construction Spending
#'   \item \code{IR} -- Import Price Index
#'   \item \code{CPILFESL} -- Core Consumer Price Index
#'   \item \code{PCEPILFE} -- Core PCE Price Index
#'   \item \code{PCEPI} -- PCE Price Index
#'   \item \code{PERMIT} -- Building Permits
#'   \item \code{TCU} -- Capacity Utilization Rate
#'   \item \code{BUSINV} -- Business Inventories
#'   \item \code{ULCNFB} -- Unit Labor Cost
#'   \item \code{IQ} -- Export Price Index
#'   \item \code{GACDISA066MSFRBNY} -- Empire State Mfg Index
#'   \item \code{GACDFSA066MSFRBPHI} -- Philadelphia Fed Mfg Index
#'   \item \code{PCEC96} -- Real Consumption Spending
#'   \item \code{GDPC1} -- Real Gross Domestic Product
#'   \item \code{ICSA} -- Weekly Unemployment Claims
#'   \item \code{DGS10} -- 10-year Treasury rates
#'   \item \code{T10Y2Y} -- 2-10 year Treasury rate spread
#'   \item \code{WALCL} -- Total Assets
#'   \item \code{PALLFNFINDEXM} -- Global Price Index of All Commodities
#'   \item \code{FEDFUNDS} -- Federal Funds Effective Rate
#'   \item \code{PPIACO} -- Producer Price Index All Commodities
#'   \item \code{CIVPART} -- Labor Force Participation Rate
#'   \item \code{M2NS} -- M2 Money Supply
#'  }
#' 
#' @return Returns the following matrices of forecasted variables:
#' \itemize{
#'  \item{\code{"interpolated_and_extrapolated"}} Returns a \code{data.frame} of the linear interpolated and \link{NNS.ARMA} extrapolated values to replace \code{NA} values in the original \code{variables} argument.  This is required for working with variables containing different frequencies, e.g. where \code{NA} would be reported for intra-quarterly data when indexed with monthly periods.
#'  \item{\code{"relevant_variables"}} Returns the relevant variables from the dimension reduction step.
#'
#'  \item{\code{"univariate"}} Returns the univariate \link{NNS.ARMA} forecasts.
#'
#'  \item{\code{"multivariate"}} Returns the multi-variate \link{NNS.reg} forecasts.
#'
#'  \item{\code{"ensemble"}} Returns the ensemble of both \code{"univariate"} and \code{"multivariate"} forecasts.
#'  }
#'
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#'
#' Viole, F. (2019) "Multi-variate Time-Series Forecasting: Nonparametric Vector Autoregression Using NNS"  \doi{10.2139/ssrn.3489550}
#'
#' Viole, F. (2020) "NOWCASTING with NNS"  \doi{10.2139/ssrn.3589816}
#'
#'
#' @examples
#'
#'  \dontrun{
#'  ## Interpolates / Extrapolates all variables to current month
#'  NNS.nowcast(h = 0)
#'  
#'  ## Additional regressors and sources specified
#'  NNS.nowcast(h = 0, additional.regressors = c("SPY", "USO"), 
#'              additional.sources = c("yahoo", "yahoo"))
#'              
#'               
#'  ### PREDICTION INTERVALS 
#'  ## Store NNS.nowcast output
#'  nns_estimates <- NNS.nowcast(h = 12)           
#'  
#'  # Create bootstrap replicates using NNS.meboot (GDP Variable)
#'  gdp_replicates <- NNS.meboot(nns_estimates$ensemble$GDPC1, 
#'                               rho = seq(0,1,.25), 
#'                               reps = 100)["replicates",]
#'                               
#'  replicates <- do.call(cbind, gdp_replicates)
#'  
#'  # Apply UPM.VaR and LPM.VaR for desired prediction interval...95 percent illustrated
#'  # Tail percentage used in first argument per {LPM.VaR} and {UPM.VaR} functions
#'  lower_GDP_CIs <- apply(replicates, 1, function(z) LPM.VaR(0.025, 0, z))
#'  upper_GDP_CIs <- apply(replicates, 1, function(z) UPM.VaR(0.025, 0, z))
#'  
#'  # View results
#'  cbind(nns_estimates$ensemble$GDPC1, lower_GDP_CIs, upper_GDP_CIs)
#'  }
#'
#' @export


NNS.nowcast <- function(h = 1,
                        additional.regressors = NULL,
                        additional.sources = NULL,
                        naive.weights = FALSE,
                        specific.regressors = NULL,
                        start.date = "2000-01-03",
                        keep.data = FALSE,
                        status = TRUE,
                        ncores = NULL){

  if(!is.null(additional.regressors) && length(additional.sources)!=length(additional.regressors)) stop("Please specify the SOURCE for each additional.regressor")
  
  variables <- c("PAYEMS", "JTSJOL",  "CPIAUCSL", "DGORDER", "RSAFS",
                 "UNRATE", "HOUST", "INDPRO", "DSPIC96", "BOPTEXP",
                 "BOPTIMP", "TTLCONS", "IR", "CPILFESL", "PCEPILFE",
                 "PCEPI", "PERMIT", "TCU", "BUSINV", "ULCNFB",
                 "IQ", "GACDISA066MSFRBNY", "GACDFSA066MSFRBPHI", "PCEC96", "GDPC1",
                 "ICSA",
                  "DGS10", "T10Y2Y", "WALCL", "PALLFNFINDEXM", "FEDFUNDS", "PPIACO", "CIVPART", "M2NS")

  sources <- c(rep("FRED", length(variables)), additional.sources)
  
  variable_list <- data.table::data.table(as.character(c(variables, additional.regressors)), sources)

  symbols <- as.character(unlist(variable_list[, 1]))
  
  if(!is.null(specific.regressors)) variable_list <- variable_list[symbols%in%symbols[specific.regressors], , drop=FALSE]

  symbols <- as.character(unlist(variable_list[, 1]))
  sources <- as.character(unlist(variable_list[, 2]))
  
  NNSdata <- new.env()
  
  for(i in 1:length(symbols)){
      quantmod::getSymbols(symbols[i], env = NNSdata, src = sources[i])
  }
  
  fetched_symbols <- ls(envir = NNSdata)
  
  if(length(fetched_symbols) < length(symbols)){
    missing_variables <- symbols[!symbols%in%fetched_symbols]
    missing_sources <- sources[!symbols%in%fetched_symbols]
    
    for(i in 1:length(missing_variables)){
        quantmod::getSymbols(missing_variables[i], src = missing_sources[i])
    }
  }
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  raw_econ_variables <- lapply(mget(symbols, envir = NNSdata), function(x) xts::to.monthly(x)[,4])
  
  
  if(!keep.data) rm(list = ls(), envir = NNSdata)
  
  econ_variables <- Reduce(function(...) merge(..., all=TRUE), raw_econ_variables)[paste0(start.date,"::")]

  colnames(econ_variables) <- symbols

  options(warn = oldw)
  
  NNS.VAR(econ_variables, h = h, tau = 12, status = status, ncores = ncores, nowcast = TRUE, naive.weights = naive.weights)
}

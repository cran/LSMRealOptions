
#' Simulate the geometric Brownian motion (GBM) stochastic process through Monte Carlo simulation
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'Simulated <- GBM.Simulate(n = 100,
#'                          t = 1,
#'                          mu = 0.05,
#'                          sigma = 0.2,
#'                          S0 = 100,
#'                          dt = 1/12)
#'
#'# ->
#'
#'## 100 simulations of 1 year of monthly price paths:
#'Simulated <- GBM_simulate(n = 100,
#'                          t = 1,
#'                          mu = 0.05,
#'                          sigma = 0.2,
#'                          S0 = 100,
#'                          dt = 1/12)
#'
#'@keywords internal
#'@export
GBM.Simulate <- function(n, t, mu, sigma, S0, dt){

#warning deprecation:
.Deprecated(msg = "'GBM.Simulate()' was deprecated in LSMRealOptions 0.2.0. \n Please use 'GBM_simulate()' instead.")

## The updated function:
return(GBM_simulate(n, t, mu, sigma, S0, dt))

}


#' Value American-Style Options Through Least-Squares Monte Carlo (LSM) Simulation
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'# Price a vanilla American put option on an asset that follows
#'# Geometric Brownian Motion
#'
#'## Step 1 - Simulate stock prices:
#'stock_prices <- GBM_simulate(n = 100, t = 1, mu = 0.05,
#'                          sigma = 0.2, S0 = 100, dt = 1/2)
#'
#'## Step 2 - Value the American put option:
#'
#'option_value <- LSM.AmericanOption(state.variables = stock_prices,
#'                                  payoff = stock_prices,
#'                                  K = 100,
#'                                  dt = 1/2,
#'                                  rf = 0.05)
#'# ->
#'option_value <- LSM_american_option(state_variables = stock_prices,
#'                                  payoff = stock_prices,
#'                                  K = 100,
#'                                  dt = 1/2,
#'                                  rf = 0.05)
#'
#'@keywords internal
#'@export
LSM.AmericanOption <- function(
  state.variables,
  payoff,
  K,
  dt,
  rf,
  call = FALSE,
  orthogonal = "Power",
  degree = 2,
  cross.product = TRUE,
  verbose = FALSE){

  #warning deprecation:
  .Deprecated(msg = "'LSM.AmericanOption()' was deprecated in LSMRealOptions 0.2.0. \n Please use 'LSM_american_option()' instead")

  return(LSM_american_option(state_variables = state.variables,
                      payoff = payoff,
                      K = K,
                      dt = dt,
                      rf = rf,
                      call = call,
                      orthogonal = orthogonal,
                      degree = degree,
                      cross_product = cross.product,
                      verbose = verbose))

}

#'Value capital investment projects through Least-Squares Monte Carlo (LSM) simulation:
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#' @examples
#'# Example: Value a capital investment project where the revenues follow a
#'# Geometric Brownian Motion stochastic process:
#'
#'## Step 1 - Simulate asset prices:
#'AssetPrices <- GBM_simulate(n = 100, t = 10, mu = 0.05,
#'                          sigma = 0.2, S0 = 100, dt = 1/2)
#'
#'## Step 2 - Perform Real Option Analysis (ROA):
#'ROA <- LSM.RealOption(state.variables = AssetPrices,
#'                      NCF = AssetPrices - 100,
#'                      CAPEX = 1000,
#'                      dt = 1/2,
#'                      rf = 0.05)
#'# ->
#'
#'## Step 2 - Perform Real Option Analysis (ROA):
#'ROA <- LSM_real_option(state_variables = AssetPrices,
#'                      NCF = AssetPrices - 100,
#'                      CAPEX = 1000,
#'                      dt = 1/2,
#'                      rf = 0.05)
#'
#'@keywords internal
#'@export
LSM.RealOption <- function(
  state.variables,
  NCF,
  CAPEX,
  dt,
  rf,
  construction = 0,
  orthogonal = "Laguerre",
  degree = 9,
  cross.product = TRUE,
  verbose = FALSE,
  debugging = FALSE){

  #warning deprecation:
  warning("'LSM.RealOption()' was deprecated in LSMRealOptions 0.2.0. \n Please use 'LSM_real_option()' instead")

  return(LSM_real_option(state_variables = state.variables,
                         NCF = NCF,
                         CAPEX = CAPEX,
                         dt = dt,
                         rf = rf,
                         construction = construction,
                         orthogonal = orthogonal,
                         degree = degree,
                         cross_product = cross.product,
                         verbose = verbose,
                         debugging = debugging))


}

#' Value operationally flexible capital investment projects through least-squares Monte Carlo (LSM) simulation:
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#' @examples
#'
#'# Example: Value a capital investment project where the revenues follow a
#'# Geometric Brownian Motion stochastic process:
#'
#'## Step 1 - Simulate asset prices:
#'asset_prices <- GBM_simulate(n = 100, t = 10, mu = 0.05,
#'                          sigma = 0.2, S0 = 100, dt = 1/2)
#'
#'## Step 2 - Perform Real Option Analysis (ROA):
#'
#'
#'ROA <- LSM.RealOption.OF(state.variables = asset_prices,
#'                      NCF = asset_prices - 100,
#'                      CAPEX = 1000,
#'                      dt = 1/2,
#'                      rf = 0.05,
#'                      Suspend.CAPEX = 100,
#'                      Suspend.OPEX = 10,
#'                      Resume.CAPEX = 100,
#'                      Abandon.CAPEX = 0
#'                      )
#'
#' # ->
#'
#'ROA <- LSM_real_option_OF(state_variables = asset_prices,
#'                      NCF = asset_prices - 100,
#'                      CAPEX = 1000,
#'                      dt = 1/2,
#'                      rf = 0.05,
#'                      suspend_CAPEX = 100,
#'                      suspend_OPEX = 10,
#'                      resume_CAPEX = 100,
#'                      abandon_CAPEX = 0
#'                      )
#'@keywords internal
#'@export
LSM.RealOption.OF <- function(
  state.variables,
  NCF,
  CAPEX,
  dt,
  rf,
  construction = 0,
  orthogonal = "Laguerre",
  degree = 9,
  cross.product = TRUE,
  Suspend.CAPEX = NULL,
  Suspend.OPEX = NULL,
  Resume.CAPEX = NULL,
  Abandon.CAPEX = NULL,
  Save.States = FALSE,
  verbose = FALSE,
  debugging = FALSE){

  #warning deprecation:
  .Deprecated(msg = "'LSM.RealOption.OF()' was deprecated in LSMRealOptions 0.2.0. \n Please use 'LSM_real_option_OF()' instead")

  ## The updated function:
  return(LSM_real_option_OF(state_variables = state.variables,
                     NCF = NCF,
                     CAPEX = CAPEX,
                     dt = dt,
                     rf = rf,
                     construction = construction,
                     orthogonal = orthogonal,
                     degree = degree,
                     cross_product = cross.product,
                     suspend_CAPEX = Suspend.CAPEX,
                     suspend_OPEX = Suspend.OPEX,
                     resume_CAPEX = Resume.CAPEX,
                     abandon_CAPEX = Abandon.CAPEX,
                     save_states = Save.States,
                     verbose = verbose,
                     debugging = debugging))
}


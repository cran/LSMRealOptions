## GBM.Simulate
#' Simulate a Geometric Brownian Motion (GBM) Stochastic Process Through Monte Carlo Simulation
#'
#' @description
#'
#' GBM is a commonly used stochastic process to simulate the price paths of stock prices and other assets, in which the log of the asset follows a random walk process with drift.
#' The \code{GBM.Simulate} function utilizes antithetic variates as a
#' simple variance reduction technique.
#'
#' @param n The total number of price paths to simulate
#' @param t The forecasting period, in years
#' @param mu The drift term of the GBM process
#' @param sigma The volatility term of the GBM process
#' @param S0 The initial value of the underlying asset
#' @param dt The discrete time step of observations, in years
#'
#'@details
#'
#'A stochastic process S(t) is a geometric brownian motion that follows the following continuous-time stochastic differential equation:
#'\deqn{ \frac{dS(t)}{S(t)} = \mu dt + \sigma dW(t)}{dS(t)/S(t) = mu dt + sigma dW(t)}
#'
#'Where \eqn{\mu}{'mu'} is the drift term, \eqn{\sigma}{'sigma'} the volatility term and \eqn{W_{t}}{W(t)} is defined as a Weiner process.
#'
#'The GBM is log-normally distributed.
#'
#'@return A matrix of simulated price paths of the GBM process. Each column corresponds to a simulated price path, and each
#'row corresponds to a simulated observed price of the simulated price paths at each discrete time period.
#'
#'@examples
#'## 100 simulations of 1 year of monthly price paths:
#'Simulated <- GBM.Simulate(n = 100,
#'                          t = 1,
#'                          mu = 0.05,
#'                          sigma = 0.2,
#'                          S0 = 100,
#'                          dt = 1/12)
#'@export
GBM.Simulate <- function(n, t, mu, sigma, S0, dt){

if(any(!is.numeric(c(n,t,mu,sigma,S0,dt)))) stop("arguments must be numeric!")
if(length(c(n,t,mu,sigma,S0,dt))>6) stop("arguments must be of length 1!")

#Calculations:
n <- round(n)
N_sim <- ifelse(n %% 2 == 0, n, n + 1)
nloops <- N_sim/2
time_periods <- seq(0, t, dt)
nsteps <- length(time_periods) - 1

##Drift:
drift <- (mu - 0.5 * sigma^2) * dt
##Shock:
shock <- matrix(stats::rnorm((nsteps) * nloops, sd = sigma), ncol = nloops) * sqrt(dt)

output <- matrix(NA, nrow = nsteps+1, ncol = N_sim)
output[1,] <- log(S0)

#Values:
output[2:(nsteps+1), seq(1, N_sim, 2)] <- log(S0) + apply(drift + shock, MARGIN = 2, cumsum)
#Antithetic Values:
output[2:(nsteps+1), seq(2, N_sim, 2)] <- log(S0) + apply(drift - shock, MARGIN = 2, cumsum)

return(exp(output[,1:n]))
}

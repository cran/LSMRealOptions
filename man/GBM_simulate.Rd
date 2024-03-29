% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDE_simulate.R
\name{GBM_simulate}
\alias{GBM_simulate}
\title{Simulate the geometric Brownian motion (GBM) stochastic process through Monte Carlo simulation}
\usage{
GBM_simulate(n, t, mu, sigma, S0, dt)
}
\arguments{
\item{n}{The total number of price paths to simulate}

\item{t}{The forecasting period, in years}

\item{mu}{The drift term of the GBM process}

\item{sigma}{The volatility term of the GBM process}

\item{S0}{The initial value of the underlying asset}

\item{dt}{The discrete time step of observations, in years}
}
\value{
A matrix of simulated price paths of the GBM process. Each column corresponds to a simulated price path, and each
row corresponds to a simulated observed price of the simulated price paths at each discrete time period.
}
\description{
GBM is a commonly used stochastic process to simulate the price paths of stock prices and other assets, in which the log of the asset follows a random walk process with drift.
The \code{GBM_simulate} function utilizes antithetic variates as a simple variance reduction technique.
}
\details{
A stochastic process S(t) is a geometric brownian motion that follows the following continuous-time stochastic differential equation:
\deqn{ \frac{dS(t)}{S(t)} = \mu dt + \sigma dW(t)}{dS(t)/S(t) = mu dt + sigma dW(t)}

Where \eqn{\mu}{'mu'} is the drift term, \eqn{\sigma}{'sigma'} the volatility term and \eqn{W_{t}}{W(t)} is defined as a Weiner process.

The GBM is log-normally distributed.
}
\examples{
## 100 simulations of 1 year of monthly price paths:
Simulated <- GBM_simulate(n = 100,
                         t = 1,
                         mu = 0.05,
                         sigma = 0.2,
                         S0 = 100,
                         dt = 1/12)
}

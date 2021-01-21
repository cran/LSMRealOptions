## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 6,
  fig.align = 'center'
)

## ----setup--------------------------------------------------------------------
library(LSMRealOptions)
set.seed(1)

## -----------------------------------------------------------------------------
# Step 1 - Simulate stock prices:
StockPrices <- GBM.Simulate(n = 1e4, t = 1, mu = 0.06, sigma = 0.2, S0 = 36, dt = 1/50)

## -----------------------------------------------------------------------------
# Step 2 - Value American put option:
PutOptionValue <- LSM.AmericanOption(state.variables = StockPrices,
                                  payoff = StockPrices,
                                  K = 40,
                                  dt = 1/50,
                                  rf = 0.06,
                                  verbose = TRUE)
print(round(unlist(PutOptionValue[1:5]),4))

## -----------------------------------------------------------------------------
# Step 1 - Simulate asset prices:
AssetPrices <- array(dim = c(51, 1e3, 2))
for(i in seq_len(2)) {
  AssetPrices[,,i] <- GBM.Simulate(n = 1e3, t = 1, mu = 0.06, 
                                  sigma = c(0.2, 0.3)[i], S0 = c(38, 35)[i], dt = 1/50)}

## -----------------------------------------------------------------------------
# Step 2 - Value American-style option:
OptionValue <- LSM.AmericanOption(state.variables = AssetPrices,
                                  payoff = pmax(AssetPrices[,,1], AssetPrices[,,2]),
                                  K = 40,
                                  dt = 1/50,
                                  rf = 0.06,
                                  verbose = TRUE,
                                  cross.product = TRUE,
                                  orthogonal = "Laguerre",
                                  degree = 9)
print(round(unlist(OptionValue[1:5]),4))

## -----------------------------------------------------------------------------
## Exercise opportunities per year:
dt <- 1/50
## strike price :
K <- 40
## short-term interest rate:
rf <- 0.06
## 100,000 simulations (50% antithetic):
N.simulations <- 1e5
## Stock price volatility:
sigma <- rep(c(rep(0.2,2),rep(0.4,2)),5)
## Stock price:
S0 <- sort(rep(seq(36,44,2),4))
## Option maturity:
TTM <- rep(1:2, 10)

LSM.output <- matrix(0, 20, 2, dimnames = list(NULL, c("Simulated American", "(s.e)")))

## Cycle through the rows of the table:
for(i in 1:20){

Simulated.values <- GBM.Simulate(n = N.simulations, t = TTM[i], 
                                 mu = rf, sigma = sigma[i], S0 = S0[i], dt = dt)

## American option pricing through LSM Simulation
output <- LSM.AmericanOption(state.variables = Simulated.values,
                   payoff = Simulated.values,
                   call = FALSE,
                   K = K,
                   dt = dt,
                   rf = rf,
                   verbose = TRUE,
                   orthogonal = "Laguerre",
                   degree = 3
                   )
LSM.output[i,1] <- output$Value
LSM.output[i,2]  <- output$`Standard Error`
}

## Compile and print results:
LnS.Table1 <- cbind.data.frame(S = S0, sigma = sigma, T = TTM, LSM.output)
print(round(LnS.Table1,3))

## -----------------------------------------------------------------------------
# Step 1 - Simulate the underlying asset price:

## Initial underlying price:
Initial.price <- 36

## discrete time step:
dt <- 1/12

## Project lifetime (in years):
Project.lifetime <- 10
forecasting.periods <- seq(0, Project.lifetime, dt)

RevenuePrices <- GBM.Simulate(n = 1e3, t = Project.lifetime, mu = 0.06, 
                              sigma = 0.2, S0 = Initial.price, dt = dt)

# Step 2 - Evaluate cash flows:

## Fixed cash flow:
FCF <- 1e4 * Initial.price

## Net cash flow is equal to variable cash flows subtract fixed cash flows:
NCF <- (1e4 * RevenuePrices - FCF) * dt

## Financial Parameters:
construction <- 0.5 / dt
rf <- 0.05

## Initial capital investment:
learning.rate <- 0.01
CAPEX <- 1e5 * exp(- learning.rate * dt * (1:nrow(RevenuePrices)-1))

# Step 3 - Evaluate Project Value through Real Options Analysis:

ProjectValue <- LSM.RealOption(state.variables = RevenuePrices,
                              NCF = NCF,
                              CAPEX = CAPEX,
                              dt = dt,
                              rf = rf,
                              construction = construction,
                              verbose = TRUE)
print(format(unlist(ProjectValue[1:6]), big.mark = ","))

## -----------------------------------------------------------------------------
## Evaluate Project Value with OF through ROA:
ProjectValue.OF <- LSM.RealOption.OF(state.variables = RevenuePrices,
                              NCF = NCF,
                              CAPEX = CAPEX,
                              dt = dt,
                              rf = rf,
                              construction = construction,
                              Suspend.CAPEX = 0.1 * CAPEX[1],
                              Suspend.OPEX = 0.05 * CAPEX[1] * dt,
                              Resume.CAPEX = 0.1 * CAPEX[1],
                              Abandon.CAPEX = 0.2 * CAPEX[1],
                              Save.States = TRUE,
                              verbose = TRUE,
                              debugging = TRUE
                              )

print(format(unlist(ProjectValue.OF[1:7]), big.mark = ","))


## -----------------------------------------------------------------------------
matplot(forecasting.periods, cbind(ProjectValue$`Cumulative Investment Prob`, 
        ProjectValue.OF$`Cumulative Investment Prob`), type = 'l', ylim = c(0,1), 
        xlab = "Forecasting Horizon", ylab = "Cumulative Investment Proportion", 
        main = "Cumulative Investment Prop. over Forecasting Horizon")
legend("right", c("ROV", "ROV + OF"),cex=0.8, col = 1:2, fill = 1:2)

## ---- warning = FALSE---------------------------------------------------------
States.list <- apply(matrix(colnames(ProjectValue.OF$`Project States`)), 1, 
                     FUN = function(x) cbind.data.frame(x, ProjectValue.OF$`Project States`[,x], 
                                                        forecasting.periods))

States.ggplot <- suppressWarnings(dplyr::bind_rows(States.list))
States.ggplot[,1] <- factor(States.ggplot[,1], levels = rev(colnames(ProjectValue.OF$`Project States`)))
colnames(States.ggplot) <- c("state", "count", "time")

library(ggplot2)

ggplot(States.ggplot, aes(x = time, y = count, fill = state)) + 
geom_bar(position = "fill", stat = "identity", width = 1) + 
scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.1)) + 
scale_x_continuous(breaks = seq(0, Project.lifetime, 1)) + 
ggtitle("Proportion of Project States over Project Planning Horizon") + xlab("Planning Horizon (Years)")

## -----------------------------------------------------------------------------
# Instantiate iterations:
it <- 0
Current.price <- Initial.price
  
# Begin Investment Trigger Value Calculate:
repeat{
  
  # Step 1: Calculate the ROV using real options analysis
  LSM_Results <- LSM.RealOption(state.variables = RevenuePrices,
                                NCF = NCF,
                                CAPEX = CAPEX,
                                dt = dt,
                                rf = rf,
                                construction = construction)

  NPV <- LSM_Results$NPV
  WOV <- LSM_Results$WOV

  # Step 2: Evaluate the next initial asset price through the 'secant' method:

  ## For the first iteration, use an arbitrary initial price multiplier of 2
  if(it == 0){
    multiplier = 2
    New.price = Current.price * multiplier
  }
  if(it > 0){

    ## NPV - a linear function of initial prices, so we can find it exactly after two iterations:
    NPV.gradient = (NPV - NPV.old) / (Current.price - Old.price)
    NPV.New.price = Current.price + (0 - NPV)/NPV.gradient
    if(it == 2) NPV.crit.value = NPV.New.price

    ## ROV -  Secant Method:
    New.price = Current.price - WOV * ((Current.price - Old.price) / (WOV - WOV.old))

    ## Which is a multiple of:
    multiplier = New.price / Current.price

    ## The WOV does not have to be exactly zero. Having it within a tolerance value 
    ## can be adequate and decrease processing time:
    WOV.tolerance <- abs(WOV) < 100
    ## If the price is identical within one cent, this can be considered the critical value:
    Price.tolerance <- round(New.price,2)==round(Current.price, 2)
    ## If the underlying asset impacts costs, and the iteration has pushed the price of the asset
    ## below zero, it's never optimal to invest immediately:
    Negative.Price <- New.price < 0
    ## Recursion break to ensure infinite loop does not occur:
    Break.loop <- it > 20
    ##Approximate the root of WOV to 2 significant figures:
    if(Price.tolerance || WOV.tolerance || Negative.Price || Break.loop){
      ROV.crit.value = New.price
      break
    } 
    }
  # Step 3: Update values:

  ## Updating simulated prices:
  RevenuePrices <- RevenuePrices * multiplier
  ## Updating the NCF of each period:
  NCF <- 1e4 * RevenuePrices - FCF
  
  ## Updating values
  Old.price <- Current.price
  Current.price <- New.price
  WOV.old <- WOV
  NPV.old <- NPV

  # Step 4: Re-iterate:
  it <- it + 1
}

print(round(c(NPV = NPV.crit.value, ROV = ROV.crit.value),2))


## -----------------------------------------------------------------------------
print(NFCP::SS.Oil$Two.Factor[2:7])

## -----------------------------------------------------------------------------
# Step 1 - List project parameters:

## Initial Price:
Initial.oil.price <- 20
## Initial State vector:
Initial.State.Vector <- c(log(Initial.oil.price), 0)

## discrete time step:
dt <- 1/12

## Project lifetime (in years):
Project.lifetime <- 10
forecasting.periods <- seq(0, Project.lifetime, dt)

## Fixed cash flow:
FCF <- 1e4 * Initial.price

## Financial Parameters:
construction <- 0.5 / dt
rf <- 0.05

CAPEX <- 1e4

#Step 1 - Simulate spot prices:
##100,000 antithetic simulations of one year of monthly observations
Simulated.Oil.Prices <- NFCP::Spot.Price.Simulate(
  X.0 = Initial.State.Vector,
  parameters = NFCP::SS.Oil$Two.Factor,
  t = 10,
  dt = dt,
  n = 1e5,
  antithetic = T,
  verbose = T)

Revenue.Prices <- Simulated.Oil.Prices$Prices

State.Variables = array(dim = c(dim(Simulated.Oil.Prices$Prices), 3))
State.Variables[,,1:2] = Simulated.Oil.Prices$State_Variables
## Include the price as a state variable:
State.Variables[,,3] = Simulated.Oil.Prices$Prices

## Net cash flow of simulated price paths:
NCF <- 1000 * Revenue.Prices - 500 * Initial.price


## -----------------------------------------------------------------------------
ProjectValue <- LSM.RealOption(state.variables = State.Variables,
                              NCF = NCF,
                              CAPEX = CAPEX,
                              dt = dt,
                              rf = rf,
                              construction = construction,
                              orthogonal = "Laguerre",
                              degree = 9,
                              verbose = T)
print(format(round(unlist(ProjectValue[1:6]),2), big.mark = ","))


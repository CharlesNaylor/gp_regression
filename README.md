# GP Regression

## Summary

This code will demonstrate the use of Gaussian processes in a dynamic linear regression, as a replacement for the Kalman Filter. More generally, Gaussian processes can be used in nonlinear regressions in which the relationship between xs and ys is assumed to vary smoothly with respect to the values of the xs. We will assume that the relationship varies smoothly with respect to time, but is static across values of xs within a given time. Another use for Gaussian processes is as a nonlinear regression technique, so that the relationship between x and y varies smoothly with respect to the values of xs, like a continuous version of random forest regressions.

You can find comments and reasoning in a set of notebooks in the /doc directory.

## Process

  1.  [Gathering Data](https://cdn.rawgit.com/billWalker/gp_regression/f34154e9/doc/Gathering_Data.html): Initial raw data retrieval and cleanup.
  2.  [Calculating Factors](https://cdn.rawgit.com/billWalker/gp_regression/9a06ddfa/doc/Calculating_Factors.html): Turning raw data into normalized factors.
  3.  [Model Overview and Simple Implementation](https://cdn.rawgit.com/billWalker/gp_regression/edac7693/doc/Specifying_the_Model-Overview_and_Simplest_Implementation.html) State the full model, then start by validating with a single factor, single asset version using known parameters.
  4.  [Model Validation - Single X, Multiple Y](https://cdn.rawgit.com/billWalker/gp_regression/edac7693/doc/Specifying_the_Model-Single_X%2C_Multiple_Y.html) Validate a single factor, multiple asset version using known parameters.
  4.  [Model Validation - Full Model](https://cdn.rawgit.com/billWalker/gp_regression/2297670c/doc/Specifying_the_Model-Full_Model.html) Validate the full multi-factor, multi-asset version using known parameters.
  5.  [Forecasting](https://rawgit.com/billWalker/gp_regression/master/doc/Forecasting.html) Test the model using the actual data.
  6.  [Backtesting] Test all historical weeks using a processor cluster
  7.  [Conclusions]

## On Currency Forecasting in General

The true model for returns in any asset class is the combination of carry and price movements. In FX, the carry is relatively stable but subject to punctuated equilibrium as new data comes to light. In developed markets, this data primarily consists of central bank rate decisions and economic data that might affect those decisions. Price movements, however, are the deterministic result of countless iterative, interacting agents. While we may know many of these agents' motives, it is not possible to aggregate their behavior with any accuracy, because the actions of each agent are affected by those of all of the other agents, and small measurement errors compound. It is impossible to tell, for example, the periodicity of market data. A day's worth of price movements at 5-minute increments looks the same as a year's worth of daily movements. The asset price measures the result of a chaotic, nonlinear dynamical system.

So why bother? What I have given you is the reason I believe no amount of layers of neural nets---or of any other algorithm which takes as an assumption that the data is stochastic---will be able to accurately predict market movements on any human timescale.  However, we know that for currencies, interest rate expectations (i.e., forecasts of carry) and other fundamental factors impact the decision-making of the agents I described above. The hope is that, over time, the small differences between predicted means, plus better-than-nothing modeling of the relationships between these factors, and the relationships between the predictions, will add up. The investor will have an edge in harvesting the carry over a portfolio of currencies, whilst weathering idiosyncratic price movements.

## Details

The Kalman filter, especially in later iterations such as the Unscented Kalman Filter or Van Der Merwe's Sigma Point Kalman filter, provides a powerful and computationally efficient method of tracking the movement of an endogenous time series given a set of correlated, but error-prone, exogenous time series. As it has a closed form solution, and operates under the Markovian assumption that D_t \perp D_{t-(2...\infty)} | D_{t-1}, i.e. that all information about the past is encapsulated in the last observation, it is particularly suitable for use in a production environment. 

However, in order to achieve this efficiency, the Kalman filter throws away anything we might learn about the nature of past intertemporal relationships, beyond what we can see in the means and covariance matrices for all variables. It also presumes Gaussian distributions. Advances on the original filter relax the Gaussian assumption, but not the Markovian one.

With a Gaussian process (GP), we can assume that parameters are related to one another in time via an arbitrary function. The disadvantage in comparison to Kalman filters is that we will wind up inverting a matrix of size T, where T is the total number of time periods in which we are interested, in order to calculate parameter values. GPs are also not necessarily solvable, and so we must rely on MCMC or its variants to evaluate the posterior distribution.

With advances in processing power, this is less of a problem than it used to be. An upcoming version of Stan is promising GPU-powered matrix inversion, and should really kick off the use of GPs in production.

As a motivating example, I'm going to use a factor regression of macroeconomic and technical indicators to forecast returns to some major currency pairs. Exploratory analysis is done in R, and the model is written in Stan. A production model would be written in Python, but I don't know that we'll get that far here.



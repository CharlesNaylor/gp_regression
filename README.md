# GP Regression

## Summary

This code will demonstrate the use of Gaussian processes in a dynamic linear regression, as a replacement for the Kalman Filter.

You can find comments and reasoning in a set of notebooks in the /doc directory.

## Process

  1.  [Gathering Data](https://cdn.rawgit.com/billWalker/gp_regression/f34154e9/doc/Gathering_Data.html): Initial raw data retrieval and cleanup.
  2.  [Calculating Factors](https://cdn.rawgit.com/billWalker/gp_regression/9a06ddfa/doc/Calculating_Factors.html): Turning raw data into normalized factors.
  3.  [Specifying the Model](https://cdn.rawgit.com/billWalker/gp_regression/9a06ddfa/doc/Specifying_the_Model-Single_Y.html)
  4.  [Model Validation]
  5.  [Forecasting]
  6.  [Model Comparison]

## Details

The Kalman filter, especially in later iterations such as the Unscented Kalman Filter or Van Der Merwe's Sigma Point Kalman filter, provides a powerful and computationally efficient method of tracking the movement of an endogenous time series given a set of correlated, but error-prone, exogenous time series. As it has a closed form solution, and operates under the Markovian assumption that D_t \perp D_{t-(2...\infty)} | D_{t-1}, i.e. that all information about the past is encapsulated in the last observation, it is particularly suitable for use in a production environment. 

However, in order to achieve this efficiency, the Kalman filter throws away anything we might learn about the nature of past intertemporal relationships, and presumes Gaussian distributions for all variables. Advances on the original filter relax the Gaussian assumption, but not the Markovian one.

With a Gaussian process (GP), we can assume that parameters are related to one another in time via an arbitrary function. The disadvantage in comparison to Kalman filters is that we will wind up inverting a matrix of size T, where T is the total number of time periods in which we are interested, in order to calculate parameter values. GPs are also not necessarily solvable, and so we must rely on MCMC or its variants to evaluate the posterior distribution.

With advances in processing power, this is less of a problem than it used to be. An upcoming version of Stan is promising GPU-powered matrix inversion, and should really kick off the use of GPs in production.

As a motivating example, I'm going to use a factor regression of macroeconomic and technical indicators to forecast returns to some major currency pairs. Exploratory analysis is done in R, the model is written in Stan, and data munging is Python. 
For the sake of reproducibility, I have tried to use only freely available data. Swap rates do not seem to be available for free, however.


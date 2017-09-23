# GP Regression

## Summary

This code will demonstrate the use of Gaussian processes in a dynamic linear regression, as a replacement for the Kalman Filter. More generally, Gaussian processes can be used in nonlinear regressions in which the relationship between xs and ys is assumed to vary smoothly with respect to the values of the xs. We will assume that the relationship varies smoothly with respect to time, but is static across values of xs within a given time. It's more usual to use Gaussian processes as a nonlinear regression technique, so that the relationship between x and y varies smoothly with respect to the values of xs, like a continuous version of random forest regressions.

You can find comments and reasoning in a set of notebooks in the /doc directory.

## Process

  1.  [Gathering Data](https://cdn.rawgit.com/billWalker/gp_regression/f34154e9/doc/Gathering_Data.html): Initial raw data retrieval and cleanup.
  2.  [Calculating Factors](https://cdn.rawgit.com/billWalker/gp_regression/9a06ddfa/doc/Calculating_Factors.html): Turning raw data into normalized factors.
  3.  [Model Overview and Simple Implementation](https://cdn.rawgit.com/billWalker/gp_regression/edac7693/doc/Specifying_the_Model-Overview_and_Simplest_Implementation.html)
  4.  [Model Validation - Single X, Multiple Y](https://cdn.rawgit.com/billWalker/gp_regression/edac7693/doc/Specifying_the_Model-Single_X%2C_Multiple_Y.html)
  4.  [Model Validation - Full Model]
  5.  [Forecasting]
  6.  [Model Comparison]

## Details

The Kalman filter, especially in later iterations such as the Unscented Kalman Filter or Van Der Merwe's Sigma Point Kalman filter, provides a powerful and computationally efficient method of tracking the movement of an endogenous time series given a set of correlated, but error-prone, exogenous time series. As it has a closed form solution, and operates under the Markovian assumption that D_t \perp D_{t-(2...\infty)} | D_{t-1}, i.e. that all information about the past is encapsulated in the last observation, it is particularly suitable for use in a production environment. 

However, in order to achieve this efficiency, the Kalman filter throws away anything we might learn about the nature of past intertemporal relationships, beyond what we can see in the means and covariance matrices for all variables. It also presumes Gaussian distributions. Advances on the original filter relax the Gaussian assumption, but not the Markovian one.

With a Gaussian process (GP), we can assume that parameters are related to one another in time via an arbitrary function. The disadvantage in comparison to Kalman filters is that we will wind up inverting a matrix of size T, where T is the total number of time periods in which we are interested, in order to calculate parameter values. GPs are also not necessarily solvable, and so we must rely on MCMC or its variants to evaluate the posterior distribution.

With advances in processing power, this is less of a problem than it used to be. An upcoming version of Stan is promising GPU-powered matrix inversion, and should really kick off the use of GPs in production.

As a motivating example, I'm going to use a factor regression of macroeconomic and technical indicators to forecast returns to some major currency pairs. Exploratory analysis is done in R, and the model is written in Stan. A production model would be written in Python, but I don't know that we'll get that far here.



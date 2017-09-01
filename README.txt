# GP Regression

This code will demonstrate the use of Gaussian processes in a dynamic linear regression, as a replacement for the Kalman Filter. The Kalman filter, especially in later iterations such as the Unscented Kalman Filter or Van Der Merwe's Sigma Point Kalman filter, provides a powerful and computationally efficient method of tracking the movement of an endogenous time series given a set of correlated, but error-prone, exogenous time series. As it has a closed form solution, and operates under the Markovian assumption that D_t \perp D_{t-(2...\infty)} | D_{t-1}, i.e. that all information about the past is encapsulated in the last observation, it is particularly suitable for use in a production environment. 

However, in order to achieve this efficiency, the Kalman filter throws away anything we might learn about the nature of past intertemporal relationships, and presumes Gaussian distributions for all variables. Advances on the original filter relax the Gaussian assumption, but not the Markovian one.

With a Gaussian process (GP), we can assume that parameters are related to one another in time via an arbitrary function. The disadvantage in comparison to Kalman filters is that we will wind up inverting a matrix of size T, where T is the total number of time periods in which we are interested, in order to calculate parameter values. GPs are also not necessarily solvable, and so we must rely on MCMC or its variants to evaluate the posterior distribution.

With advances in processing power, this is less of a problem than it used to be. An upcoming version of Stan is promising GPU-powered matrix inversion, and should really kick off the use of GPs in production.

As a motivating example, I'm going to use a factor regression of macroeconomic and technical indicators to forecast returns to some major currency pairs. Exploratory analysis is done in R, the model is written in Stan, and data munging is Python. 
For the sake of reproducibility, I will try to use only freely available data.

You can find comments and reasoning in a set of notebooks in the /doc directory.

Regression-based-Backtesting

The following code was written for my bachelors thesis on "Regression-based Backtesting of Volatility Forecasts".

In this work we evaluate and compare forecast optimality tests for MZ regressions via Backtesting when performing volatility forecasts with parametric forecast models.
Our focus in this work was on the classic Mincer-Zarnowitz (MZ) regression and its extensions, which involves the transformation of the dependent and independent variables.

In the MZ regression, the variable of interest - in our case, volatility - is regressed against a constant and the forecasted value. After conducting these regressions, various MZ tests for forecast optimality can be applied. In our work, we considered three different tests:

1. Classic MZ test

2. Global F-test (which does not consider the intercept)

3. F-test that includes the intercept.

We conducted a Monte Carlo simulation to evaluate the size and power of these three tests. To do this, we generated data following a GARCH model and used this data to perform forecasts from a predefined forecasting model in which we can change the persistence. Subsequently, we applied MZ regressions and the respective MZ tests. In our simulation, we repeated this procedure 1000 times to compute the size and power of the tests.

Since volatility is latent, we used volatility proxies in our simulation. Besides using squared returns as proxies, we also computed realized variances.
We created a function that allows for artificial generation of realized variances with different numbers of intra daily returns.

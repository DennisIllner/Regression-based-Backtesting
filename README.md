# Regression-based-Backtesting

#The following code was written for my bachelors thesis on "Regression-based Backtesting of Volatility Forecasts".

#In this work we evaluate and compare forecast optimality tests for MZ regressions via Backtesting when performing volatility forecasts with parametric forecast models
#Our focus in this work was on the classic Mincer-Zarnowitz (MZ) regression and its extensions, which involve the transformation of both dependent and independent variables.

#In the MZ regression, the variable of interest - in our case, volatility - is regressed against a constant and the forecasted value. After conducting these regressions, #various MZ tests for forecast optimality can be applied. In our work, we considered three different tests:

1. Classic MZ test

2. Global F-test (which does not consider the intercept)

3. F-test that includes the intercept.

#ferre

FERRE is a data analysis code written in FORTRAN90. It matches models to data, taking a
set of observations and identifying the model parameters that best reproduce the data, in a chi-
squared sense. Model predictions are to be given as an array whose values are a function of the
model parameters, i.e. numerically. FERRE holds this array in memory, or in a direct-access
binary file, and interpolates in it. The code returns, in addition to the optimal set of parameters,
their error covariance, and the corresponding model prediction.

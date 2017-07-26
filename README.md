# ODonnelletal_2017_imbalances
Model code for O'Donnell et al, 2017.
Study on multidimensional imbalances in neural circuit properties in brain disorders.

Files include:

BRIAN CIRCUIT MODEL OF L2/3 SOMATOSENSORY CORTEX (PYTHON)
- layer23model.py; Brian simulation that sets up single L2/3 circuit model and simulates one pattern of input activity, plots summary results.
- layer23model_varyparams.py; Brian simulation of same circuit model that averages over many runs, input patterns, and input sparsities, then varies each parameter by 20%.

LOGISTIC MODEL (MATLAB)
- dfr_dcorr_dtheta.m; Numerically compute the gradients in firing rate and correlation with respect to the threshold and slope of the logistic model.
- errfunc.m; Computes error function, needed for gradient descent in model fitting.
- fr_corr_logistic.m; Returns expected firing rate and correlation from pair of identical neurons with logistic input-output functions parameterized by threshold and slope parameters.
- get_alpha_logistic.m; Returns offset parameter of logistic function, given the slope and x-y location of any point on the curve.
- get_basicstats_prob_3param.m; Calculates the mean firing rate, s.d. in firing rates, and mean pairwise correlation in a raster movie of firing rates. Assumes a Poisson firing model.
- mcmc_logistic_3param.m; Fits the five parameters of a population logistic model to three activity statistics of the data, using a markov-chain monte carlo method for gradient descent.
- samples_logistic.m; Generates binary sample data from 5-parameter population logistic model.
- stats_from_logistic_params.m; Calculates the 3 target stats (firing rate mean and s.d., pairwise covariance mean) from logistic model parameters. Involves doing numerical integration over implied probability distributions.
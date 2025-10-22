# Changepoint-Scoring

This repo contains all my code related to my project on "Change-point Scoring". 
Broadly, how can we apply score estimation to improve test for change-points and 
change-point estimates. 

The repo is structured into the following folders:
- "R": contains all helper functions, not to be used by the user including 
common functions related to noise distributions, calculation of change-point
estimators, and functions run by simulations
- "simulations": contains the main "run_simulations.R" script which runs a borad
range of simulations with progress / parallel processing and config files used
to run simulations
- "results": saved output of simulation runs
- "data": saved plots / historical simulation runs (old code)
- "analysis": contains functions to analyse output of simulation runs / other
simple analysis files

Work in progress.

-------------------------------------------------------------------------------

Ideas (theory)
 - change-point scoring as a transform of the data by rho_hat()/Fisher_information,
 theoretical results about standard CUSUM statistics can then be applied based on
 this transformed data
    - okay when transform function estimated independently from the data (eg: 
    via sample splitting), proving results may be tricky when score estimated
    from all the data
- Is the above transform the "best" (in what way) transform onto Gaussian data.
LAN links?? Is mean 0 and variance 1

Ideas (code)
- What if we learn the score of the data with a CP initially rather than trying 
to estimate CP and delta, with a fixed initial score?
- Implement new delta estimation via scoring difference / Fisher information
- Check if it is actually Fisher Information (ie: rho' applied to eps not Y) [true
either way for Gaussian data / t-dist]

-------------------------------------------------------------------------------
The Scoring Algorithm:
Initialise: data, Y | initial (normalised) score, $\rho$ | 

1. Estimate CP via tau_hat =  max(scaling * CUSUM(rho(Y)))
2. Estimate delta_hat =
mean_post_tau(rho(Y)) - mean_pre_tau(rho(Y)) (recursively at minimal cost)
3. Compute residuals as deCPed Y (not necessarily mean 0)
4. Update rho (normalised ie: divide by an estimate of the FI if necessary)

Q: Can we estimate the normalised rho()/I(rho) directly / fit directly to the data?
Does ASM do this?

-------------------------------------------------------------------------------


$$
\hat{\delta} = \frac{\sum_{t= \tau + 1}^{n} \hat{\rho}(Y_t) - \sum_{t= 1}^{\tau} \hat{\rho}(Y_t)}{I(\rho)}
$$
  
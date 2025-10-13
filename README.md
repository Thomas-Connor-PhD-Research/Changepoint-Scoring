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
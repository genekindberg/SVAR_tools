# SVAR_tools
Structural VAR (SVAR) toolbox for bayesian VAR estimation and a range of identification methods

VAR toolbox that allows a bayesian estimation of the VAR coefficients with changeable priors (edit in BVAR in "Functions" folder). 

Allows multiple identifications: 
* Spectral and Limited Spectral restrictions - identifying the shock that maximizes the share of variance in a desired frequency band - see [Dieppe, Francis, and Kindberg-Hanlon 2021 (JEDC)](https://www.sciencedirect.com/science/article/pii/S0165188921001512) and [Dieppe, Francis, and Kindberg-Hanlon 2021 (ECB Working Paper)](https://www.ecb.europa.eu/pub/pdf/scpwps/ecb.wp2534~2383e60ba4.en.pdf?eeffb1db5c1033da86491dbc6c96ce9a)
* The long-run restriction of Blanchard & Quah 1989/Gali 1999
* The Max Share approach of Francis et. al. 2014
* A Cholesky identification.
* Sign and zero restrictions
* Restrictions on the share of FEVD explained by the different shocks.


The bayesian VAR function outputs:
* Impulse responses (IRFs)
* Forecast error variance decompositions (FEVDs)
* Historical decompositions
* Structural shock series

Full examples are provided.

## Long-run and Max Share restrictions
The file *RunMain.m* estimates a VAR consisting of US labor productivity, employment, investment and consumption as a share of GDP, inflation, and long-term bond yields. It then demonstrates how to estimate a technology shock using the Spectral, Limited Spectral, long-run, and Max Share restrictions and plots them side by side.


## Sign, zero, and relative FEVD contributions
The file *RunMain_signandzero.m* estimates a VAR consisting of US labor productivity, employment, inflation, and long-term bond yields. It then demonstrates a method of identifying a technology shock, demand shock, monetary policy shock, and a supply shock using sign and zero restrictions. It also demonstrates how to impose FEVD magnitude restrictions. For example, technology shocks are assumed to have a larger share of the variance of labor productivity at the 5 year horizon, while demand shocks are assumed to dominate the FEVD of labor productivity in the first year. This is for demonstrative purposes rather than for any particular theoretical reason.



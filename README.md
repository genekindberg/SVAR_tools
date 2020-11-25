# SVAR_tools
Structural VAR (SVAR) toolbox for bayesian VAR estimation and a range of identification methods

VAR toolbox that allows a bayesian estimation of the VAR coefficients with changeable priors (edit in BVAR in "Functions" folder). 

Allows multiple identifications: 
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
The file *RunMain.m* estimates a VAR consisting of US labor productivity, employment, investment and consumption as a share of GDP, inflation, and long-term bond yields. It then demonstrates how to use long-run and Max Share restrictions to identify technology shocks and plots them side by side.

## Sign, zero, and relative FEVD contributions
The file *RunMain_signandzero.m* estimates a VAR consisting of US labor productivity, employment, inflation, and long-term bond yields. It then demonstrates a method of identifying a technology shock, demand shock, monetary policy shock, and a supply shock using sign and zero restrictions. It also demonstrates how to impose FEVD magnitude restrictions. For example, technology shocks are assumed to have a larger share of the variance of labor productivity at the 5 year horizon, while demand shocks are assumed to dominate the FEVD of labor productivity in the first year. This is for demonstrative purposes rather than for any particular theoretical reason.

Further identifications will be provided in due course, including spectral and limited spectral identifications. Futher information available in: "The Identification of Dominant Macroeconomic Drivers: Coping with Confounding Shocks", available at genekindberg.github.io/Research 

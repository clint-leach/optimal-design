

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5172768.svg)](https://doi.org/10.5281/zenodo.5172768)



# optimal-design

Example code for using recursive Bayesian computation in optimal design, using simulated data and binary regression. 

## Publication

Leach, Clinton B., Perry J. Williams, Joseph M. Eisaguirre, Jamie N. Womble, Michael R. Bower, and Mevin B. Hooten, submitted 2021, Recursive Bayesian computation facilitates adaptive optimal design in ecological studies, *Ecology*

## File list

`output/data.rds`  
`scripts/plots.R`  
`src/designs.jl`  
`src/example.jl`  
`src/probit.reg.jl`  
`optimal-design.Rproj`  
`LICENSE`  
`README.md`

## Description

`output/data.rds` – R data file containing a data frame of simulated data with columns `loc` – integer site index; `y` – binary presence/absence at the site; `x` – covariate value at the site; `K1` – binary indicator of inclusion in the initial sample; `score` – computed design criteria for including that site in a second data collection effort

`scripts/plots.R` – R script for generating Figure 1 of the above manuscript

`src/designs.jl` – Julia code defining a function to loop over designs and predicted data, apply PPRB, and compute quantities needed to evaluate the design criteria

`src/example.jl` – Julia script simulating binary data, fitting an initial binary regression model, and computing the design criteria for each of the potential designs

`src/probi.reg.jl` – Julia code defining functions to fit an initial probit regression model to a complete data set (probit_reg) and to use PPRB to update an initial model fit given new data (probit_pprb)

`optimal-design.Rpro`j – Rproject file for loading the project in Rstudio

`LICENSE` – file specifying the license of the included code

`README.md` – file containing same metadata

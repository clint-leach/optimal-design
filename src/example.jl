using Distributions
using Plots
using PDMats
using BenchmarkTools
using LinearAlgebra
using StatsPlots
using DataFrames
using RCall

include("probit.reg.jl")
include("designs.jl")

# Size of full domain
N = 100

# Size of initial sample
n = 10

# Parameters
β = [0.0, 0.5]

# Data
X = hcat(fill(1.0, N), sort(rand(Normal(0, 3), N)))
p = cdf.(Normal(), X * β)
y = rand.(Bernoulli.(p))

# Selecting observations
obs = sample(1:N, n, replace = false)

# Fitting model
chain = probit_reg(y[obs], X[obs, :], X, 500000)

# Optimal design ==========================================================

# Multiple imputation
ppred = mean(chain.y, dims = 2)
ypred = hcat(fill(0, N), fill(1, N))

# Designs
D = findall(x -> !(x in obs), 1:N)
L = length(D)

# PPRB with post-hoc MI
scores_pprb = designs_pprb(ypred, D, chain.p, chain.β)

vars = sum(scores_pprb.vars, dims = 1)[1, :, :]
meanvar = ppred[D] .* vars[:, 2] + (1.0 .- ppred[D]) .* vars[:, 1]

# Output data for plotting in R ================================================

data = DataFrame(loc = 1:100, y = y, x = X[:, 2], K1 = 0, score = NaN)
data[obs, :K1] = 1
data[D, :score] = meanvar

R"
saveRDS($(chain.p), 'output/chain.rds')
saveRDS($(data), 'output/data.rds')
"

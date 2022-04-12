using DifferentialEquations
using DiffEqPhysics
using PlotlyJS
using DataFrames
import Statistics as stat
using StaticArrays
using LoopVectorization
using LsqFit
using BenchmarkTools
using Distributions

# Import Constants
include("utils/constants.jl");

# Import Functions
include("utils/get_ystart_v4.jl")
include("utils/hsimd.jl")
include("utils/store_orbit.jl")
include("utils/calcxyzae.jl")
include("utils/gettt.jl")
include("utils/calcmeanperiods.jl")
include("utils/plotting_functions.jl")
include("utils/readtt.jl")
include("utils/calc_nbody.jl")
include("utils/likelihood.jl")
include("utils/logprior.jl")

# Import Planet System parameters
include("utils/KOI2433.jl");
# include("utils/LHS1678.jl") 

sol = calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan);

#Store model
@time df = store_orbit(sol, mass, nbody, NVEC);
#Extract TTs
@time nTT,TT=getTT(sol, mass, G);

#Calcalate Mean Periods
tt_period,tt_T0=calcmeanperiods(nTT,TT);

#Do a comparision of the period we asked for vs the period we got.
println(tt_period ./ DAYS)
println(periods ./ DAYS)
println((tt_period ./ DAYS) - (periods ./ DAYS))

# makeplots_v3(df,nbody,planet_names)

plotTTVs(tt_T0, tt_period, nTT, TT)

modelpars = [mass; periods; T0; sqecosω; sqesinω];
extrapars = [ep, tspan, G, nbody];
ll = likelihood(modelpars, extrapars)

priors = [mass_prior; periods_prior; T0_prior; sqecosω_prior;  sqesinω_prior]
lp = logprior(modelpars, priors)


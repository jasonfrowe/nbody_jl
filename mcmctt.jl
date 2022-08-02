using DifferentialEquations
using DiffEqPhysics
#using PlotlyJS
using Plots
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
include("utils/likelihood_v2.jl")
include("utils/logprior_v2.jl")
include("utils/dmcmc.jl")
include("utils/betarescale.jl")
include("utils/plot_1ttv.jl")

# Import Planet System parameters
include("utils/KOI2433.jl");
# include("utils/LHS1678.jl") 
# include("utils/V1298tau.jl")

nTT_obs
TT_obs
TT_errobs

Δt=0.00012
sol=calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)
data=store_orbit(sol, mass::AbstractVector, nbody::Integer, NVEC::Integer)
transit_data=getTT(sol, mass::AbstractVector, G::AbstractFloat)
nTT=transit_data[1]
TT=transit_data[2]
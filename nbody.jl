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
include("utils/dmcmc.jl")
include("utils/betarescale.jl")

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

# makeplots_v3(df,nbody,planet_names) #Generate N-body plot.
 
# Plot TTVs
plotTTVs(tt_T0, tt_period, nTT, TT)

modelpars = [mass; periods; T0; sqecosω; sqesinω];
extrapars = [ep, tspan, G, nbody];
ll = likelihood(modelpars, extrapars)

# Set up Matrix wit Priors 
priors = [mass_prior; periods_prior; T0_prior; sqecosω_prior;  sqesinω_prior];
lp = logprior(modelpars, priors)

# Set up β parameter for Gibbs sampling
β = [massβ; periodsβ; T0β; sqecosωβ; sqesinωβ]; 

llx = ll + lp
llx1, x1, ac = mhgmcmc(modelpars, extrapars, llx, β, priors);

nsample = 1000
# TODO@jasonfrowe add threads to genchain for multiple walkers
chain = genchain(modelpars, extrapars, β, priors, nsample);

# get likelihood of last chain to check that the walker is working.
likelihood(chain[nsample+1][1], extrapars)

burnin = 100
corscale = betarescale(chain[nsample+1][1], extrapars, β, nsample, burnin)

nsample = 5000
chain = genchain(modelpars, extrapars, β .* corscale, priors, nsample);

#Calculuate acceptance rates

nfitpars=length(modelpars)
accrate=zeros(nfitpars)
accsel=zeros(nfitpars)
vecrate=0
vecsel=0
for i in 2:nsample
    if chain[i][2][2] > 0
        if chain[i][2][1]==0
            accrate[Int(chain[i][2][2])]+=1
        end
        accsel[Int(chain[i][2][2])]+=1
    else
        if chain[i][2][1]==0
            vecrate+=1
        end
        vecsel+=1
    end 
end
println(accrate./accsel)
println(accsel)
print([vecrate/vecsel,vecrate,vecsel])



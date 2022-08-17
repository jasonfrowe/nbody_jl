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
using LinearAlgebra
using ProgressMeter

1+1

# Import Constants
include("utils/constants.jl");

# Import Functions
include("utils/get_ystart_v4.jl")
include("utils/get_ystart_v5.jl")
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
include("utils/plot_ttv.jl")

# Import Planet System parameters
include("utils/KOI2433.jl");
# include("utils/LHS1678.jl") 
# include("utils/V1298tau.jl")
# include("utils/testcase.jl")

#Get model solution
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

modelpars = [mass; periods; T0; sqecosω; sqesinω];
extrapars = [ep, tspan, G, nbody];
ll = likelihood(modelpars, extrapars)

nb = 2
plot_1ttv(nTT_obs, TT_obs, TTerr_obs, nTT, TT, modelpars, extrapars, nb)

# makeplots_v3(df,nbody,planet_names) #Generate N-body plot.
 
# Plot TTVs
plotTTVs(tt_T0, tt_period, nTT, TT)

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

plot_TTV(chain, modelpars, extrapars, TT_obs, TTerr_obs, nTT_obs, 2)

# get likelihood of last chain to check that the walker is working.
likelihood(chain[nsample+1][1], extrapars) .+ logprior(chain[nsample+1][1], priors)

burnin = 100
corscale = betarescale(chain[nsample+1][1], extrapars, β, nsample, burnin)

modelpars_cstart = chain[nsample+1][1] .* 1.0

likelihood(modelpars_cstart, extrapars) .+ logprior(modelpars_cstart, priors)

nsample = 50000
chain = genchain(modelpars_cstart, extrapars, β .* corscale, priors, nsample);

chain2 =  genchain(modelpars_cstart, extrapars, β .* corscale, priors, nsample, nsample, chain)

nsample3 = 500000
chain3 =  genchain(chain2[nsample+1][1], extrapars, β .* corscale, priors, nsample3, nsample, chain2)


likelihood(chain3[nsample3+1][1], extrapars) .+ logprior(chain3[nsample3+1][1], priors)

#Calculuate acceptance rates

nfitpars=length(modelpars)
accrate=zeros(nfitpars)
accsel=zeros(nfitpars)
vecrate=0
vecsel=0
for i in 2:nsample
    if chain3[i][2][2] > 0
        if chain3[i][2][1]==0
            accrate[Int(chain3[i][2][2])]+=1
        end
        accsel[Int(chain3[i][2][2])]+=1
    else
        if chain3[i][2][1]==0
            vecrate+=1
        end
        vecsel+=1
    end 
end
println(accrate./accsel)
println(accsel)
print([vecrate/vecsel,vecrate,vecsel])

plot_TTV(chain3, modelpars, extrapars, TT_obs, TTerr_obs, nTT_obs, 2)
plot_TTV(chain3, modelpars, extrapars, TT_obs, TTerr_obs, nTT_obs, 3)
plot_TTV(chain3, modelpars, extrapars, TT_obs, TTerr_obs, nTT_obs, 4)
plot_TTV(chain3, modelpars, extrapars, TT_obs, TTerr_obs, nTT_obs, 5)

maximum(TTerr_obs)
indexin(maximum(TTerr_obs),TTerr_obs)

mplot = zeros(nsample)
for i ∈ 1:nsample
    mplot[i] = chain3[i][1][2] / MEARTH
end


trace = histogram(x = mplot, nbinsx=50);
plot(trace)


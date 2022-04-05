using DifferentialEquations
using DiffEqPhysics
using PlotlyJS
using DataFrames
import Statistics as stat
using StaticArrays
using LoopVectorization
using LsqFit

# Import Constants
include("utils/constants.jl");

# Import Functions
include("utils/get_ystart.jl")
include("utils/get_ystart_v2.jl")
include("utils/hsimd.jl")
include("utils/store_orbit.jl")
include("utils/calcxyzae.jl")
include("utils/gettt.jl")
include("utils/calcmeanperiods.jl")
include("utils/plotting_functions.jl")

# Import Planet System parameters
include("utils/KOI2433.jl");

tspan = (tstart, tend); #integration range
p0,q0=get_ystart_v2(mass,eccn,periods,T0,ep); #getting initial conditions

q0=SVector{size(q0)[1]}(q0); #store as Static Array for better speed
p0=SVector{size(p0)[1]}(p0);
mass=SVector{size(mass)[1]}(mass);
param=mass;

prob = HamiltonianProblem(H_simd, p0, q0, tspan, param); #set up Hamiltonian to solve

#Calculate Model
Δt=0.00012; #minimum(periods[2:nbody])/10
@time sol_001 = solve(prob, KahanLi8(), abstol=1e-10, reltol=1e-10, dt = Δt);

#Store model
@time df_001 = store_orbit(sol_001, mass, nbody, NVEC);

#Extract TTs
@time nTT,TT=getTT(sol_001, mass, G);

#Calcalate Mean Periods
tt_period,tt_T0=calcmeanperiods(nTT,TT);

println(tt_period ./ DAYS)
println(periods ./ DAYS)
#println(tt_T0 ./ DAYS)
println((tt_period ./ DAYS) - (periods ./ DAYS))


#makeplots_v3(df_001,nbody,planet_names)

plotTTVs(tt_T0,tt_period,nTT,TT)

t=0
using BenchmarkTools
@benchmark H_simd(p0, q0, mass, t)


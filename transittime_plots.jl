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

##################################################The first nested for loop is to calculate the cnages in transit times for all planets in the system#################################
########################The second for loop is just finding the limit of planet f since its transit times fluctuate more drastically as Δt changes####################################

#Initial values for the dataframes based on value of Δt=0.00012
Δt=0.00012
sol=calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)
data=store_orbit(sol, mass::AbstractVector, nbody::Integer, NVEC::Integer)
transit_data=getTT(sol, mass::AbstractVector, G::AbstractFloat)
nTT=transit_data[1]
TT=transit_data[2]

nTT_len=size(nTT)[1]


p2=plot()


Δt_ar2=[Δt]
for planet ∈ 2:nTT_len
    p2=plot()
    raw_data=[]
    for i ∈ 1:10
        new_Δt=Δt+0.000005
        Δt=new_Δt
        push!(Δt_ar2, Δt)

        #Updating dataframes given the Δt change
        sol=calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)
        transit_data=getTT(sol, mass::AbstractVector, G::AbstractFloat)
        nTT=transit_data[1]
        TT=transit_data[2]./DAYS .*86400 
    
        nTT_len=size(nTT)[1]

        planet_data=nTT[planet]
        data_point=TT[planet, planet_data]
        push!(raw_data, data_point)
    end
    fixed_data=raw_data .- mean(raw_data)
    
    for i ∈ 1:10
        plot!(p2,[Δt_ar2[i]], [fixed_data[i]], seriestype = :scatter, legend = false, title = "Transit times for varying (+)Δt for planet "*string(planet_names[planet]), minorgrod=true)
        xlabel!("Value of Δt")
        ylabel!("Transit Times (seconds)") 
    end
    display(p2)
end


#Finding Δt limit for planet f
#5E-17 is the limit, any smaller than that doesn't plot properly. 
#Changes can be seen in transit time on order of 1E-5 seconds
#Changed in Δt of 5E-10 cause changed in the transit time on order of 5-10 seconds

Δt=0.00012

pf=plot()
raw_data=[]
Δt_arf=[Δt]
for i ∈ 1:25
    new_Δt=Δt+5E-10
    Δt=new_Δt
    push!(Δt_arf, Δt)

    #Updating dataframes given the Δt change
    sol=calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)
    transit_data=getTT(sol, mass::AbstractVector, G::AbstractFloat)
    nTT=transit_data[1]
    TT=transit_data[2]./DAYS .*86400
    
    nTT_len=size(nTT)[1]

    planet_data=nTT[6]
    data_point=TT[6, planet_data]
    push!(raw_data, data_point)
end
fixed_data=raw_data .- mean(raw_data)
    
for i ∈ 1:10
    plot!(pf,[Δt_arf[i]], [fixed_data[i]], seriestype = :scatter, legend = false, title = "Transit times for varying (+5E-10)Δt for planet f ", minorgrid=true)
    xlabel!("Value of Δt")
    ylabel!("Transit Times (seconds)") 
end
display(pf)

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


#Initial values for the dataframes based on arbitrary value of Δt=0.00012
Δt=0.00012
sol=calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)
data=store_orbit(sol, mass::AbstractVector, nbody::Integer, NVEC::Integer)

TotalEnergy = data.TotalEnergy
sizeTE=size(TotalEnergy)[1]
firstTE=TotalEnergy[1]
lastTE=TotalEnergy[sizeTE]
difference=abs(firstTE)-abs(lastTE)

p2=plot()
plot!(p2, [Δt], [difference], marker = :rect, color = :red, seriestype = :scatter, legend = false, title = "Change in Total Energy for varying (+)Δt", minorgrid=true)
xlabel!("Value of Δt")
ylabel!("Change in Energy")

#Not necessary to explore 0<Δt<0.00012 because we found that the changes were very negligible and don't start diverging until t>0.00012

#Loop for Δt>0.00012
diff_ar2=[]
Δt_ar2=[Δt]
for i ∈ 1:27
    new_Δt=Δt+0.000005
    Δt=new_Δt
    push!(Δt_ar2, Δt)

    #Updating dataframes given the Δt change
    sol=calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)
    data=store_orbit(sol, mass::AbstractVector, nbody::Integer, NVEC::Integer)
    
    #Extracting relevant arrays for calculating
    TotalEnergy = data.TotalEnergy
    
    #Finding the change in total energy and appending it to an array to store the values
    sizeTE=size(TotalEnergy)[1]
    firstTE=TotalEnergy[1]
    lastTE=TotalEnergy[sizeTE]
    difference=abs(firstTE)-abs(lastTE)
    push!(diff_ar2, difference)

    #Plotting instructions
    plot!(p2, [Δt_ar2[i]], [diff_ar2[i]], seriestype = :scatter, legend = false, title = "Change in Total Energy for varying (+)Δt", minorgrid=true)
    xlabel!("Value of Δt")
    ylabel!("Change in Energy")
end
display(p2)


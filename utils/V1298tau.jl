planet_names=['A', 'b', 'c', 'd', 'e'];

#mrexp = 1.0 / 0.42 #mass Radius relation

#Particle Masses
mass = [1.101 * MSUN, 
        48.0 * MEARTH, 
        1 * MEARTH, 
        1 * MEARTH, 
        44 * MEARTH
       ]

mass_prior = [
      mass[1] 0.5 * MSUN 2.0 * MSUN 0.02 * MSUN 0 ;
      mass[2] 0.1 * MEARTH  500.0 * MEARTH 50.0 * MEARTH 3 ;
      mass[3] 0.1 * MEARTH  100.0 * MEARTH 0 1 ; 
      mass[4] 0.1 * MEARTH  100.0 * MEARTH 0 1 ;
      mass[5] 0.1 * MEARTH  500.0 * MEARTH 60.0 * MEARTH 3 
]

mβ = 1.0 * MEARTH
massβ = mβ .* ones(length(mass))

#Orbital Eccentricity
sqecosω = [0.0, 0.0, 0.0, 0.0, 0.0] #eccentricity

sqecosω_prior = [
      sqecosω[1] -1.0 1.0 0 0;
      sqecosω[2] -1.0 1.0 0 1;
      sqecosω[3] -1.0 1.0 0 1;
      sqecosω[4] -1.0 1.0 0 1;
      sqecosω[5] -1.0 1.0 0 1
]

ecβ = 0.001
sqecosωβ = ecβ .* ones(length(mass))

sqesinω = [0.0, 0.0, 0.0, 0.0, 0.0]

sqesinω_prior = [
      sqesinω[1] -1.0 1.0 0 0;
      sqesinω[2] -1.0 1.0 0 1;
      sqesinω[3] -1.0 1.0 0 1;
      sqesinω[4] -1.0 1.0 0 1;
      sqesinω[5] -1.0 1.0 0 1
]

esβ = 0.001
sqesinωβ = esβ .* ones(length(mass))

#Orbital Periods
periods = [0.0, 
           24.14042 * DAYS,
           8.24872 * DAYS, 
           12.40214 * DAYS, 
           43.36682 * DAYS
          ]  

pσ = 1.0e-4
periods_prior = [
      periods[1] 0.0 1.0 0 0;
      periods[2] 0.0 0.0 pσ 2;
      periods[3] 0.0 0.0 pσ 2;
      periods[4] 0.0 0.0 pσ 2;
      periods[5] 0.0 0.0 pσ 2
]

pβ = 1.0e-5
periodsβ = pβ .* ones(length(mass))

nbody=length(mass) #number of bodies

#Centre of transit times
T0 = [0.0, 
      8298.2095 * DAYS, 
      8293.341 * DAYS, 
      8287.804 * DAYS, 
      8310.892 * DAYS
     ]

T0σ = 3.0e-3
T0_prior = [
      T0[1] 0.0 1.0 0 0;
      T0[2] 0.0 1.0 T0σ 2;
      T0[3] 0.0 1.0 T0σ 2;
      T0[4] 0.0 1.0 T0σ 2;
      T0[5] 0.0 1.0 T0σ 2
]

Tβ = 1.0e-5
T0β = Tβ .* ones(length(mass))

#Time span for simulation
tstart = 7060.0 * DAYS #starting time
tend   = 9550.0 * DAYS #end time

#Integration time-step
Δt=0.00012; #minimum(periods[2:nbody])/10

#Read in transit times -- this uses the DR25 Kepler Database format
# ttdir = "/Volumes/astro/Kepler/Kepler_n/timing/"
ttdir = "TTVs/"
ttbase = "V1298tau"
nTT_obs, TT_obs, TTerr_obs = readtt(nbody,ttdir,ttbase);

# Things you probably do not want to change
ep = tstart #Epoch -- if you need to shift T0 to match tstart
nbody = length(mass); #compute the number of bodies 
tspan = (tstart, tend); #integration range

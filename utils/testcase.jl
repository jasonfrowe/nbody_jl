planet_names=['A', 'b', 'c', 'd'];

mrexp = 1.0 / 0.42 #mass Radius relation

#Particle Masses
mass = [1.00 * MSUN, 
        1.0 * MEARTH, 
        1.0 * MEARTH, 
        1.0 * MEARTH
       ]

mass_prior = [
      mass[1] 0.5 * MSUN 2.0 * MSUN 0.02 * MSUN 0 ;
      mass[2] 0.1 * MEARTH  50.0 * MEARTH 0 1 ;
      mass[3] 0.1 * MEARTH  50.0 * MEARTH 0 1 ; 
      mass[4] 0.1 * MEARTH  50.0 * MEARTH 0 1 
]

mβ = 1.0 * MEARTH
massβ = mβ .* ones(length(mass))

#Orbital Eccentricity
sqecosω = [0.0, 0.0, 0.0, 0.0] #eccentricity

sqecosω_prior = [
      sqecosω[1] -1.0 1.0 0 0;
      sqecosω[2] -1.0 1.0 0 1;
      sqecosω[3] -1.0 1.0 0 1;
      sqecosω[4] -1.0 1.0 0 1
]

ecβ = 0.001
sqecosωβ = ecβ .* ones(length(mass))

sqesinω = [0.0, 0.0, 0.0, 0.0] #eccentricity

sqesinω_prior = [
      sqesinω[1] -1.0 1.0 0 0;
      sqesinω[2] -1.0 1.0 0 1;
      sqesinω[3] -1.0 1.0 0 1;
      sqesinω[4] -1.0 1.0 0 1
]

esβ = 0.001
sqesinωβ = esβ .* ones(length(mass))

#Orbital Periods
periods = [0.0, 
           44.85436722412415 * DAYS,
           85.13889843078842 * DAYS, 
           130.20806895403814 * DAYS
          ]  

pσ = 1.0e-5
periods_prior = [
      periods[1] 0.0 1.0 0 0;
      periods[2] 0.0 0.0 pσ 2;
      periods[3] 0.0 0.0 pσ 2;
      periods[4] 0.0 0.0 pσ 2
]

pβ = 1.0e-5
periodsβ = pβ .* ones(length(mass))

nbody=length(mass) #number of bodies

#Centre of transit times
T0 = [0.0, 
      814.5915799108891 * DAYS, 
      825.5277398760087 * DAYS, 
      795.9202994695297 * DAYS
     ]

T0σ = 1.0e-4
T0_prior = [
      T0[1] 0.0 1.0 0 0;
      T0[2] 0.0 1.0 T0σ 2;
      T0[3] 0.0 1.0 T0σ 2;
      T0[4] 0.0 1.0 T0σ 2
]

Tβ = 1.0e-5
T0β = Tβ .* ones(length(mass))

#Time span for simulation
tstart =  285.3977955338*DAYS #starting time
tend   = 1524.0011050734*DAYS #end time

#Integration time-step
Δt=0.00012; #minimum(periods[2:nbody])/10

#Read in transit times -- this uses the DR25 Kepler Database format
# ttdir = "/Volumes/astro/Kepler/Kepler_n/timing/"
ttdir = "TTVs/"
ttbase = "koi2433"
#nTT_obs, TT_obs, TTerr_obs = readtt(nbody,ttdir,ttbase);

# Things you probably do not want to change
ep = 780.0 #Epoch -- if you need to shift T0 to match tstart
nbody = length(mass); #compute the number of bodies 
tspan = (tstart, tend); #integration range

planet_names=['A', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'];

mrexp = 1.0 / 0.42 #mass Radius relation

#Particle Masses
mass = [0.99 * MSUN, 
        2.405^mrexp * MEARTH, 
        2.312^mrexp * MEARTH, 
        2.431^mrexp * MEARTH, 
        1.918^mrexp * MEARTH, 
        1.204^mrexp * MEARTH, 
        1.445^mrexp * MEARTH, 
        2.246^mrexp * MEARTH, 
        1.201^mrexp * MEARTH
       ]

mass_prior = [
      mass[1] 0.5 * MSUN 2.0 * MSUN 0.02 * MSUN 0 ;
      mass[2] 0.1 * MEARTH  50.0 * MEARTH 0 1 ;
      mass[3] 0.1 * MEARTH  50.0 * MEARTH 0 1 ; 
      mass[4] 0.1 * MEARTH  50.0 * MEARTH 0 1 ;
      mass[5] 0.1 * MEARTH  50.0 * MEARTH 0 1 ;
      mass[6] 0.1 * MEARTH  50.0 * MEARTH 0 1 ;
      mass[7] 0.1 * MEARTH  50.0 * MEARTH 0 1 ;
      mass[8] 0.1 * MEARTH  50.0 * MEARTH 0 1 ;
      mass[9] 0.1 * MEARTH  50.0 * MEARTH 0 1
]

mβ = 1.0 * MEARTH
massβ = mβ .* ones(length(mass))

#Orbital Eccentricity
sqecosω = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #eccentricity

sqecosω_prior = [
      sqecosω[1] -1.0 1.0 0 0;
      sqecosω[2] -1.0 1.0 0 1;
      sqecosω[3] -1.0 1.0 0 1;
      sqecosω[4] -1.0 1.0 0 1;
      sqecosω[5] -1.0 1.0 0 1;
      sqecosω[6] -1.0 1.0 0 1;
      sqecosω[7] -1.0 1.0 0 1;
      sqecosω[8] -1.0 1.0 0 1;
      sqecosω[9] -1.0 1.0 0 1
]

ecβ = 0.001
sqecosωβ = ecβ .* ones(length(mass))

sqesinω = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #eccentricity

sqesinω_prior = [
      sqesinω[1] -1.0 1.0 0 0;
      sqesinω[2] -1.0 1.0 0 1;
      sqesinω[3] -1.0 1.0 0 1;
      sqesinω[4] -1.0 1.0 0 1;
      sqesinω[5] -1.0 1.0 0 1;
      sqesinω[6] -1.0 1.0 0 1;
      sqesinω[7] -1.0 1.0 0 1;
      sqesinω[8] -1.0 1.0 0 1;
      sqesinω[9] -1.0 1.0 0 1
]

esβ = 0.001
sqesinωβ = esβ .* ones(length(mass))

#Orbital Periods
periods = [0.0, 
           1.5162232294E+01 * DAYS,
           1.0043770976E+01 * DAYS, 
           5.6414651546E+01 * DAYS, 
           2.7904178405E+01 * DAYS, 
           6.4674455089E-01 * DAYS, 
           6.0632603880E+00 * DAYS, 
           8.6431613712E+01 * DAYS, 
           3.3737863179E+00 * DAYS
          ]  

pσ = 1.0e-5
periods_prior = [
      periods[1] 0.0 1.0 0 0;
      periods[2] 0.0 0.0 pσ 2;
      periods[3] 0.0 0.0 pσ 2;
      periods[4] 0.0 0.0 pσ 2;
      periods[5] 0.0 0.0 pσ 2;
      periods[6] 0.0 0.0 pσ 2;
      periods[7] 0.0 0.0 pσ 2;
      periods[8] 0.0 0.0 pσ 2;
      periods[9] 0.0 0.0 pσ 2
]

pβ = 1.0e-5
periodsβ = pβ .* ones(length(mass))

nbody=length(mass) #number of bodies

#Centre of transit times
T0 = [0.0, 
      7.7470258617E+01 * DAYS, 
      2.8661761746E+02 * DAYS, 
      9.6129584036E+01 * DAYS, 
      8.7939989352E+01 * DAYS, 
      6.5026398489E+01 * DAYS, 
      6.9179060589E+01 * DAYS, 
      3.5603662008E+02 * DAYS, 
      2.8834623274E+02 * DAYS
     ]

T0σ = 1.0e-4
T0_prior = [
      T0[1] 0.0 1.0 0 0;
      T0[2] 0.0 1.0 T0σ 2;
      T0[3] 0.0 1.0 T0σ 2;
      T0[4] 0.0 1.0 T0σ 2;
      T0[5] 0.0 1.0 T0σ 2;
      T0[6] 0.0 1.0 T0σ 2;
      T0[7] 0.0 1.0 T0σ 2;
      T0[8] 0.0 1.0 T0σ 2;
      T0[9] 0.0 1.0 T0σ 2
]

Tβ = 1.0e-5
T0β = Tβ .* ones(length(mass))

#Time span for simulation
tstart =  285.3977955338*DAYS #starting time
tend   = 1524.0011050734*DAYS #end time

#Integration time-step
Δt=0.0012; #minimum(periods[2:nbody])/10

#Read in transit times -- this uses the DR25 Kepler Database format
# ttdir = "/Volumes/astro/Kepler/Kepler_n/timing/"
ttdir = "TTVs/"
ttbase = "koi2433"
nTT_obs, TT_obs, TTerr_obs = readtt(nbody,ttdir,ttbase);

# Things you probably do not want to change
ep = tstart #Epoch -- if you need to shift T0 to match tstart
nbody = length(mass); #compute the number of bodies 
tspan = (tstart, tend); #integration range

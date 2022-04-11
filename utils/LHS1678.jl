planet_names=['A', 'b', 'c', 'd'];

mrexp = 1.0 / 0.42 #mass Radius relation

#Particle Masses
mass = [0.345 * MSUN, 
        0.696^mrexp * MEARTH, 
        0.982^mrexp * MEARTH,
        0.931^mrexp * MEARTH 
       ]

#Orbital Eccentricity
sqecosω = [0.0, 0.0, 0.0, 0.0] #eccentricity
sqesinω = [0.0, 0.0, 0.0, 0.0] #eccentricity

#Orbital Periods
periods = [0.0, 
           0.8602322000E+00 * DAYS,
           3.6942470000E+00 * DAYS, 
           4.9652070000E+00 * DAYS
          ]  

nbody=length(mass) #number of bodies

#Centre of transit times
T0 = [0.0, 
      1.4114768054E+03 * DAYS, 
      1.4147594110E+03 * DAYS, 
      1.4145629000E+03 * DAYS
     ]

#Time span for simulation
tstart = 1410.0000000000*DAYS #starting time
tend   = 2122.0000000000*DAYS #end time

#Integration time-step
Δt=0.000012; #minimum(periods[2:nbody])/10

#Read in transit times -- this uses the DR25 Kepler Database format
ttdir = "/Volumes/astro/Kepler/Kepler_n/timing/"
ttbase = "koi2433"
nTT_obs, TT_obs, TTerr_obs = readtt(nbody,ttdir,ttbase);

# Things you probably do not want to change
ep = tstart #Epoch -- if you need to shift T0 to match tstart
nbody = length(mass); #compute the number of bodies 
tspan = (tstart, tend); #integration range

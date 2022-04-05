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

#Orbital Eccentricity
sqecosω = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #eccentricity
sqesinω = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #eccentricity

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

tstart =  285.3977955338*DAYS #starting time
tend   = 1524.0011050734*DAYS #end time

ep = tstart
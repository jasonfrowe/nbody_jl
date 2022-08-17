#Threads
#Threads 
nthreads=Threads.nthreads() #Record number of threads.  If nbody>10*nthreads, it might be worth using multithreading.

#Physical Constants
const GR_MKS     = 6.67408e-11 #Graviational constant SI 
const RJUP_MKS   = 71490000.0 #Radius of Jupiter (m)
const RSUN_MKS   = 695508000.0 #Radius of the Sun (m)
const MJUP_MKS   = 1.266865349218008E+08/GR_MKS*10^9 #mass of Jupiter kg
const MSUN_MKS   = 1.3271244004193938E+11/GR_MKS*10^9 #Mass of the Sun kg
const MEARTH_MKS = 3.9860043543609598E+05/GR_MKS*10^9 #Mass of the Earth kg
const AU_MKS     = 149597870700.0 #Astronomical Unit m
const DAYS_MKS   = 86400.0 #length of day (s)
const YR_MKS     = 31557600.0e0 #use 1-year = 365.25 
const π2         = π*π # π²     

# Scaled Units
const G  = 1.0    #Graviational Constant
const MU = MSUN_MKS  #mass units
const LU = AU_MKS  #Length Scale 
const TU = sqrt(LU * LU * LU / (MU * GR_MKS)) #Time scale

# #Scaled Units
# const G  =  0.00029591220819207774   #Graviational Constant
# const MU = MSUN_MKS  #mass units  
# const TU = DAYS_MKS    #time units
# const LU = (MU * TU * TU * GR_MKS)^(1.0 / 3.0) #Length scale

#Collision/Softening Scale (not currently used)
const RCOL = 2 * RJUP_MKS / LU
const SOFTENING = 0.0

#Simulation parameters
const NVEC = 6 #number of parameters to describe the state of the particle

#Mass of Earth relative to the Sun
const MSUN   = MSUN_MKS / MU
const MEARTH = MEARTH_MKS / MU
const MJUP   = MJUP_MKS / MU
const DAYS   = DAYS_MKS / TU
const YR     = YR_MKS / TU 
const AU     = AU_MKS / LU
const RSUN   = RSUN_MKS / LU
const RJUP   = RJUP_MKS / LU

#Threads 
# nthreads=Threads.nthreads() #Record number of threads.  If nbody>10*nthreads, it might be worth using multithreading.

#Physical Constants
# const RJUP_MKS   = 71490000.0 #Radius of Jupiter
# const MJUP_MKS   = 1.898e27 #mass of Jupiter kg
# #const MSUN_MKS   = 1.98911e30 #mass of Sun kg
# const MSUN_MKS = 1.9884754159665356e+30
#const YR_MKS     = 31557600.0e0 #length of year in sec
# const YR_MKS = 31558149.7635
# const GR_MKS     = 6.67408e-11 #Graviational constant
# const MEARTH_MKS = 5.97219e24 #Mass of the Earth kg
# const π2         = π*π

# #Scaled Units
# const G  = 1.0    #Graviational Constant
# const MU = MSUN_MKS  #mass units  
# const TU = YR_MKS    #time units
# const LU = (MU * TU * TU * GR_MKS)^(1.0 / 3.0) #Length scale

# #Collision/Softening Scale (not currently used)
# const RCOL = 2 * RJUP_MKS / LU
# const SOFTENING = 0.0

# #Simulation parameters
# const NVEC = 6 #number of parameters to describe the state of the particle

# #Mass of Earth relative to the Sun
# const MSUN   = MSUN_MKS / MU
# const MEARTH = MEARTH_MKS / MU
# const MJUP   = MJUP_MKS / MU
# const DAYS   = 86400.0 / TU
# const YR     = YR_MKS / TU 

#This defines 1AU
# const AU_MKS = ( (1 + MEARTH / MSUN) / (4 * π2))^(1 / 3) * LU
# const AU = AU_MKS / LU;

# #Threads 
# nthreads=Threads.nthreads() #Record number of threads.  If nbody>10*nthreads, it might be worth using multithreading.

# #Physical Constants
# RJUP_MKS   = 71490000.0 #Radius of Jupiter
# MJUP_MKS   = 1.898e27 #mass of Jupiter kg
# MSUN_MKS   = 1.98911e30 #mass of Sun kg
# YR_MKS     = 31557600.0e0 #length of year in sec
# GR_MKS     = 6.67408e-11 #Graviational constant
# MEARTH_MKS = 5.97219e24 #Mass of the Earth kg
# π2         = π*π

# #Scaled Units
# G  = 1.0    #Graviational Constant
# MU = MSUN_MKS  #mass units  
# TU = YR_MKS    #time units
# LU = (MU * TU * TU * GR_MKS)^(1.0 / 3.0) #Length scale

# #Collision/Softening Scale (not currently used)
# RCOL = 2 * RJUP_MKS / LU
# SOFTENING = 0.0

# #Simulation parameters
# NVEC = 6 #number of parameters to describe the state of the particle

# #Mass of Earth relative to the Sun
# MSUN   = MSUN_MKS / MU
# MEARTH = MEARTH_MKS / MU
# MJUP   = MJUP_MKS / MU
# DAYS   = 86400.0 / TU
# YR     = YR_MKS / TU 

# #This defines 1AU
# AU_MKS = ( (1 + MEARTH / MSUN) / (4 * π2))^(1 / 3) * LU
# AU = AU_MKS / LU;
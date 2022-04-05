#Threads 
nthreads=Threads.nthreads() #Record number of threads.  If nbody>10*nthreads, it might be worth using multithreading.

#Physical Constants
const RJUP_MKS   = 71490000.0 #Radius of Jupiter
const MJUP_MKS   = 1.898e27 #mass of Jupiter kg
const MSUN_MKS   = 1.98911e30 #mass of Sun kg
const YR_MKS     = 31557600.0e0 #length of year in sec
const GR_MKS     = 6.67408e-11 #Graviational constant
const MEARTH_MKS = 5.97219e24 #Mass of the Earth kg
const π2         = π*π

#Scaled Units
const G  = 1.0    #Graviational Constant
const MU = MSUN_MKS  #mass units  
const TU = YR_MKS    #time units
const LU = (MU * TU * TU * GR_MKS)^(1.0 / 3.0) #Length scale

#Collision/Softening Scale (not currently used)
const RCOL = 2 * RJUP_MKS / LU
const SOFTENING = 0.0

#Simulation parameters
const NVEC = 6 #number of parameters to describe the state of the particle

#Mass of Earth relative to the Sun
const MSUN   = MSUN_MKS / MU
const MEARTH = MEARTH_MKS / MU
const MJUP   = MJUP_MKS / MU
const DAYS   = 86400.0 / TU
const YR     = YR_MKS / TU 

#This defines 1AU
const AU_MKS = ( (1 + MEARTH / MSUN) / (4 * π2))^(1 / 3) * LU
const AU = AU_MKS / LU;
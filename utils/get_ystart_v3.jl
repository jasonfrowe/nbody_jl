function get_ystart_v3(mass, eccn, periods, T0, ep)
    #sets up the initial conditions of the simulation translating m,e,p,t0 -> x,y,z,vx,vy,vz
   
    nbody = size(mass)[1]
    nvecd2 = Integer(NVEC / 2)
    q = zeros(nbody * nvecd2)
    p = zeros(nbody * nvecd2)

    trueanom = zeros(nbody)
    for i in 1:nbody
        if periods[i] > 0
            ph = ((ep - T0[i]) / periods[i] - floor((T0[i] - ep)/periods[i])) * 2 * π
            trueanom[i] = ph
        else
            trueanom[i] = 0
        end
    end

    tmass = 0.0
    isort = sortperm(periods) #iterate starting for inner most object
    
    #Track centre of mass
    com = [0.0, 0.0, 0.0]
    #com_old = [0.0, 0.0, 0.0]
    
    for i ∈ isort
        asemi = (periods[i] * periods[i] * G * (mass[i] + tmass) / (4 * π2))^(1.0 / 3.0)
    
        r = asemi * (1 - eccn[i] * eccn[i]) / (1 + eccn[i] * cos(trueanom[i]))   

        x = r * cos(trueanom[i]) + com[1] #X-position
        y = r * sin(trueanom[i]) + com[2]#Y-position
        z = 0.0 + com[3] #ignoring inclination for the moment.
        
        #Update CoM
        com[1] = (x * mass[i] + com[1] * tmass) / (mass[i] + tmass)
        com[2] = (y * mass[i] + com[2] * tmass) / (mass[i] + tmass) 
        com[3] = (z * mass[i] + com[3] * tmass) / (mass[i] + tmass) 

        if periods[i] <= 0 #central object will likely not have an initial period
            temp = 0.0
        else
            temp = 2 * π * asemi / (periods[i] * sqrt(1 - eccn[i]^2))
        end
        vx = -temp * sin(trueanom[i]) #X velocity
        vy = +temp * (eccn[i] + cos(trueanom[i])) #Y velocity
        vz = 0.0

        q[i*nvecd2-2] = x
        q[i*nvecd2-1] = y
        q[i*nvecd2  ] = z
        
        p[i*nvecd2-2] = vx * mass[i]
        p[i*nvecd2-1] = vy * mass[i]
        p[i*nvecd2  ] = vz * mass[i]

        tmass += mass[i]
    end
    
    return p,q
end;
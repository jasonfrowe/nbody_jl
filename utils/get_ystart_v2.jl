function get_ystart_v2(mass,eccn,periods,T0,ep)
    #sets up the initial conditions of the simulation translating m,e,p,t0 -> x,y,z,vx,vy,vz
   
    nbody=size(mass)[1]
    nvecd2=Integer(NVEC/2)
    q=zeros(nbody*nvecd2)
    p=zeros(nbody*nvecd2)

    trueanom=zeros(nbody)
    for i in 1:nbody
        if periods[i]>0
            ph=((ep-T0[i])/periods[i]-floor((T0[i]-ep)/periods[i]))*2*pi
            trueanom[i]=ph
        else
            trueanom[i]=0
        end
    end

    #asemi = (periods.*periods*G.*(mass.+mass[1])/(4*pi2)).^(1.0/3.0)
    tmass=0.0
    isort = sortperm(periods)
    asemi=zeros(nbody)
    for i in isort
        asemi[i]=(periods[i]*periods[i]*G*(mass[i]+tmass)/(4*π2)).^(1.0/3.0)
        tmass+=mass[i]
    end
    
    r = asemi.*(1 .- eccn.*eccn)./(1 .+ eccn.*cos.(trueanom))
    
    bodies=[]
    tmass=0.0
    for i in 1:nbody
        μmass = tmass / (tmass + mass[i])
        μmass_star = mass[i] / (tmass + mass[i])
        x = r[i]*cos(trueanom[i]) #X-position
        y = r[i]*sin(trueanom[i]) #Y-position
        z = 0.0
        
        if periods[i] <= 0 #central object will likely not have an initial period
            temp = 0.0
        else
            temp = 2*π*asemi[i]/(periods[i]*sqrt(1-eccn[i]^2))
        end
        vx = - temp*sin(trueanom[i]) #X velocity
        vy = + temp*(eccn[i]+cos(trueanom[i])) #Y velocity
        vz=0.0
        
        q[i*nvecd2-2] = x*μmass
        q[i*nvecd2-1] = y*μmass
        q[i*nvecd2  ] = z*μmass
        
        p[i*nvecd2-2] = vx*mass[i]*μmass
        p[i*nvecd2-1] = vy*mass[i]*μmass
        p[i*nvecd2  ] = vz*mass[i]*μmass

        #update star

        q[1] += -x*μmass_star
        q[2] += -y*μmass_star
        q[3] += -z*μmass_star

        p[1] +=  -vx*mass[1]*μmass_star
        p[2] +=  -vy*mass[1]*μmass_star
        p[3] +=  -vz*mass[1]*μmass_star

        tmass+=mass[i]
        
    end
    
    return p,q
end;
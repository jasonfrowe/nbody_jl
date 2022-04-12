function get_ystart_v4(
        mass::AbstractVector, 
        periods::AbstractVector, 
        T0::AbstractVector, 
        ep::AbstractFloat, 
        sqecosω::AbstractVector, 
        sqesinω::AbstractVector)
    #sets up the initial conditions of the simulation translating m,e,p,t0 -> x,y,z,vx,vy,vz
   
    nbody = length(mass)
    nvecd2 = Integer(NVEC / 2)
    q = zeros(nbody * nvecd2) #position vector
    p = zeros(nbody * nvecd2) #momentum vector

    isort = sortperm(periods) #iterate starting with inner most object
    #println("isort: $(isort)")

    ϵₛ = eps(1.0) #A small number

    pid2 = π/2.0
    tpi  = 2.0*π 
    fourpisq = 4.0*π*π
    fourpid3 = 4.0*π/3.0

    eccn=(sqecosω .* sqecosω .+ sqesinω .* sqesinω) #eccentricity
    ecosω=sqrt.(eccn) .* sqecosω
    esinω=sqrt.(eccn) .* sqesinω

    ω = zeros(nbody)
    asemi = zeros(nbody)
    f = zeros(nbody)
    Omrad = zeros(nbody)
    mtot = zeros(nbody) #cumulative mass
    irad = ones(nbody) .* pid2 #inclination -- for future use
    for ii ∈ 1:nbody #loop over bodies in order of orbital period
        i=isort[ii]
        
        if eccn[i] > 1.0
            eccn[i] = 0.99 #set upper limit for useful eccentricity
        end

        if abs(eccn[i]) < ϵₛ
            ω[i] = 0.0
        else
            if abs(ecosω[i]) < ϵₛ
                ω[i] = pid2
            else
                ω[i] = atan(esinω[i]/ecosω[i])
            end

            if ecosω[i] > 0.0 && esinω[i] < 0.0
                ω[i] = tpi + ω[i]
            elseif ecosω[i] < 0.0 && esinω[i] >= 0.0
                ω[i] = π + ω[i]
            elseif ecosω[i] <= 0.0 && esinω[i] < 0.0
                ω[i] = π + ω[i]
            end
        end

        #mtot += mass[i]
        if i == 1 #star is always first as P=0
            mtot[i]=mass[i]
        else
            mtot[i]=mass[i]+mtot[isort[ii-1]]
        end

        if i == 1
            asemi[i] = 0.0
            f[i] = 0.0

        else

            asemi[i] = (G * mtot[i] * periods[i] * periods[i] / (fourpisq))^(1.0 / 3.0)
            varpi = Omrad[i] + ω[i]

            bige = 2.0 * atan(sqrt((1.0 - eccn[i]) / (1.0 + eccn[i])) * tan(0.5 * (pid2 - ω[i]))) #Murray and Dermott eqn 2.46
            Meanlong = varpi + bige - eccn[i] * sin(bige) # eqn 2.52
            Meananom = Meanlong - varpi                   # M&D 2.53

            eccanom = Meananom + tpi * (ep - T0[i]) / periods[i]
            lame = eccanom + varpi

            while lame < -1.0 * π
                lame += tpi
            end
            while lame > 1.0 * π
                lame -= tpi
            end

            bigm = lame - varpi
            bige = bigm

            de = 1.0e30 # initialize with a big number
            maxiter = 250
            iter =0 

            while((de > ϵₛ) && (iter < maxiter)) # iterate until we reach numerical precision
                bigen = bigm + eccn[i] * sin(bige) #Murray and Dermott eqn 2.54, iterative solution
                de = abs(bigen - bige)
                bige = bigen
                iter += 1
            end

            # for j ∈ 1:25
            #     bige = bigm+eccn[i] * sin(bige) 
            # end

            f[i] = 2.0 * atan(sqrt((1.0 + eccn[i]) / (1.0 - eccn[i])) * tan(0.5 * bige))
    
        end
    
    end
    # println("mtot: $(mtot)")

    r=zeros(nbody) 
    x0=zeros(nbody) ;y0=zeros(nbody); z0=zeros(nbody)
    vx0=zeros(nbody) ; vy0=zeros(nbody) ; vz0=zeros(nbody)
    x1=zeros(nbody) ; y1=zeros(nbody) ; z1=zeros(nbody)
    vx1=zeros(nbody) ; vy1=zeros(nbody) ; vz1=zeros(nbody)
    x2=zeros(nbody) ; y2=zeros(nbody) ; z2=zeros(nbody)
    vx2=zeros(nbody) ; vy2=zeros(nbody) ; vz2=zeros(nbody)
    x3=zeros(nbody) ; y3=zeros(nbody) ; z3=zeros(nbody)
    vx3=zeros(nbody) ; vy3=zeros(nbody) ; vz3=zeros(nbody)

    for ii ∈ 1:nbody #loop over bodies in order of orbital period
        i=isort[ii]

        if i > 1

            r[i] = asemi[i] * (1.0 - eccn[i] * eccn[i]) / (1.0 + eccn[i] * cos(f[i])) #M+D eqn 2.20
            x0[i] = r[i] * cos(f[i])
            y0[i] = r[i] * sin(f[i])
            z0[i] = 0.0

            vx0[i] = -sqrt(G * mtot[i] / (asemi[i] * (1.0 - eccn[i] * eccn[i]))) * sin(f[i]) # 2.36
            vy0[i] =  sqrt(G * mtot[i] / (asemi[i] * (1.0 - eccn[i] * eccn[i]))) * (eccn[i] + cos(f[i]))
            vz0[i] = 0.0

            # rotate in plane due to omega;
            x1[i] = x0[i] * cos(ω[i]) - sin(ω[i]) * y0[i] #P1, eqn 2.119
            y1[i] = x0[i] * sin(ω[i]) + cos(ω[i]) * y0[i]
            z1[i] = z0[i]

            vx1[i] = vx0[i] * cos(ω[i]) - sin(ω[i]) * vy0[i]
            vy1[i] = vx0[i] * sin(ω[i]) + cos(ω[i]) * vy0[i]
            vz1[i] = vz0[i]

            # rotate out of sky-plane due to inclination; P2, eqn 2.119
            x2[i] = x1[i] #the long axis of the ellipse stays the same
            y2[i] = y1[i] * cos(irad[i]) - z1[i] * sin(irad[i]) #the y extent of the ellipse shrinks with rotation.
            z2[i] = y1[i] * sin(irad[i]) + z1[i] * cos(irad[i]) #the z extent should grow with rotation

            vx2[i] = vx1[i]
            vy2[i] = vy1[i] * cos(irad[i]) - vz1[i] * sin(irad[i])
            vz2[i] = vy1[i] * sin(irad[i]) + vz1[i] * cos(irad[i])

            #'Precess' plane around reference axes, P3, eqn 2.121, make negative as in line 74 setup.pro
            x3[i] = -1.0*((x2[i]*cos(Omrad[i]) - y2[i]*sin(Omrad[i])))
            y3[i] = -1.0*(x2[i]*sin(Omrad[i]) + y2[i]*cos(Omrad[i]))
            z3[i] = -1.0*(z2[i])

            vx3[i] = -1.0*(vx2[i]*cos(Omrad[i]) - vy2[i]*sin(Omrad[i]))
            vy3[i] = -1.0*(vx2[i]*sin(Omrad[i]) + vy2[i]*cos(Omrad[i]))
            vz3[i] = -1.0*(vz2[i])

        end

    end

    x4=zeros(nbody)  ; y4=zeros(nbody)  ; z4=zeros(nbody)
    vx4=zeros(nbody) ; vy4=zeros(nbody) ; vz4=zeros(nbody)

    for ii ∈ 1:nbody
        i=isort[ii]
        
        x4[i]=x3[i]
        y4[i]=y3[i]
        z4[i]=z3[i]
        vx4[i]=vx3[i]
        vy4[i]=vy3[i]
        vz4[i]=vz3[i]

        for jj ∈ i:ii-1
            
            j=isort[jj]

            x4[i]=x4[i]+x3[j]*mass[j]/mtot[j]
            y4[i]=y4[i]+y3[j]*mass[j]/mtot[j]
            z4[i]=z4[i]+z3[j]*mass[j]/mtot[j]
            vx4[i]=vx4[i]+vx3[j]*mass[j]/mtot[j]
            vy4[i]=vy4[i]+vy3[j]*mass[j]/mtot[j]
            vz4[i]=vz4[i]+vz3[j]*mass[j]/mtot[j]
        end

        q[i*nvecd2-1] = x4[i] #position vector
        q[i*nvecd2-0] = -y4[i]
        q[i*nvecd2-2] = -z4[i]

        p[i*nvecd2-1] = vx4[i] * mass[i] #momentum vector
        p[i*nvecd2-0] = -vy4[i] * mass[i]
        p[i*nvecd2-2] = -vz4[i] * mass[i]

    end

    return p,q

end
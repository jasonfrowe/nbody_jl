function get_ystart_v5(
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
    #fourpisq = 4.0*π*π
    #fourpid3 = 4.0*π/3.0

    eccn=(sqecosω .* sqecosω .+ sqesinω .* sqesinω) #eccentricity
    ecosω=sqrt.(eccn) .* sqecosω
    esinω=sqrt.(eccn) .* sqesinω

    ω = zeros(nbody)
    asemi = zeros(nbody)
    f = zeros(nbody)
    Omrad = zeros(nbody)
    mtot = zeros(nbody) #cumulative mass - replaced with mu 
    irad = ones(nbody) .* pid2 #inclination -- for future use
     
    n = zeros(nbody) # 2 π / period

    mtot[1] = mass[1] #assume that first body is the central star.
    for ii ∈ 2:nbody
        i = isort[ii]
        im1 = isort[ii-1]
        mtot[i] = mtot[im1] + mass[i] #cumulative mass including current planet
    end

    for ii ∈ 2:nbody
        i = isort[ii]
        n[i] = tpi / periods[i] #periods
        n2 = n[i]*n[i]
        asemi[i] = (G*mtot[i]/n2)^(1.0/3.0) #semi-major axis
    end

    for ii ∈ 2:nbody
        i = isort[ii]
        if eccn[i] > ϵₛ
            ω[i] = atan(esinω[i]/eccn[i],ecosω[i]/eccn[i]) 
        end
        if ω[i] < 0.0
            ω[i] += tpi
        end
    end

    MA = zeros(nbody)
    for ii ∈ 2:nbody
        i = isort[ii]
        theta = π/2.0
        varpi = ω[i] + Omrad[i] #never used! 
        f1 = theta - ω[i]  #TrueAnomaly
        tanE_over_2 = tan(f1/2.0)*sqrt((1.0-eccn[i])/(1.0+eccn[i]))
        E1 = 2.0*atan(tanE_over_2)
        MeanAnomaly = E1 - eccn[i]*sin(E1)
        MA[i] = MeanAnomaly - n[i] * (T0[i] - ep)
        if MA[i] < 0
            MA[i] += tpi
        end
    end

    E = zeros(nbody)
    convcrit = 1.0e-15
    stateraw = zeros(nbody,NVEC)
    for ii ∈ 2:nbody
        i = isort[ii]

        E[i] = MA[i]*1.0
        δE = 1
        while abs(δE) > convcrit
            δE = (E[i] - eccn[i] * sin(E[i]) - MA[i]) / (1.0 - eccn[i] * cos(E[i]))
            E[i] = E[i] - δE
        end

        tanf_over_2 = tan(E[i]/2.0)*sqrt( (1.0+eccn[i]) / (1.0-eccn[i]) )
        f[i] = 2.0 * atan(tanf_over_2)
        r = asemi[i] * (1.0 - eccn[i] * cos(E[i]))
        x = r * cos(f[i])
        y = r * sin(f[i])
        xdot = -n[i] * asemi[i] * sin(f[i]) / sqrt( 1.0 - eccn[i]*eccn[i])
        ydot = n[i] * asemi[i] * (eccn[i] + cos(f[i])) / sqrt(1.0 - eccn[i]*eccn[i])

        P1 = [ cos(ω[i]) -sin(ω[i]) 0; 
            sin(ω[i]) cos(ω[i]) 0; 
            0 0 1 ]
        P2 = [ 1 0 0; 
            0 cos(irad[i]) -sin(irad[i]);
            0 sin(irad[i]) cos(irad[i]) ]
        P3 = [ cos(Omrad[i]) -sin(Omrad[i]) 0; 
            sin(Omrad[i]) cos(Omrad[i]) 0; 
            0 0 1 ]

        xyz = [x; y; 0]
        vxyz = [xdot; ydot; 0]

        #Matrix multiplication
        Y1 = similar(xyz) ; mul!(Y1, P1, xyz)
        Y2 = similar(Y1) ; mul!(Y2, P2, Y1)
        XYZ = similar(Y2) ; mul!(XYZ, P3, Y2)
        VY1 = similar(vxyz) ; mul!(VY1, P1,vxyz)
        VY2 = similar(VY1) ; mul!(VY2, P2, VY1)
        VXYZ = similar(VY2) ; mul!(VXYZ, P3, VY2)

        stateraw[i,:] = [ XYZ[1] XYZ[2] XYZ[3] VXYZ[1] VXYZ[2] VXYZ[3] ]

    end

    state = zeros(nbody,NVEC)
    for ii ∈ 2:nbody
        i = isort[ii]
        state[i,:] = stateraw[i,:] .* 1.0
        for jj ∈ 1:ii-1
            j = isort[jj]
            #println([i,j])
            state[i,:] += stateraw[j,:] * mass[j] / (mtot[j]) 
        end
        #println(' ')
    end

    for i ∈ 1:nbody
        q[i*nvecd2-2] =  state[i,3]
        q[i*nvecd2-1] = -state[i,1] 
        q[i*nvecd2-0] = -state[i,2]

        p[i*nvecd2-2] =  state[i,6] * mass[i]
        p[i*nvecd2-1] = -state[i,4] * mass[i]
        p[i*nvecd2-0] = -state[i,5] * mass[i]
    end


    return p,q

end
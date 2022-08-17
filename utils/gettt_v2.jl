function getTT_v2(sol, mass::AbstractVector, G::AbstractFloat)
    #Calculates transit-times for the system.
    #assumes mass[1]=central star.  Transit occurs when y=0 and x>0.


    nbody = length(mass) #number of bodies
    nvecd2 = Integer(NVEC / 2) #dimensions (3)
    nvd2bd = nvecd2 * nbody 
    #nvbd = nvec * nbody
    npt = size(sol.u)[1] #number of samples

    TT = zeros(nbody,npt) #store the transit times
    nTT = zeros(Int64,nbody) #the number of transit times for each planet

    SR = zeros(nvecd2) #Star position
    SV = zeros(nvecd2) #Star Velocity
    p = zeros(nvd2bd) #momentum
    q = zeros(nvd2bd) #position
    # Rvec = zeros(nvecd2) #current position vector
    # Rvec_old = zeros(nbody,nvecd2) #store reprevious position vector
    t_old = 0.0 #store previous time stamp

    rstar = 1.0 * RSUN #Placeholder for size of star.

    tcmax = 5 #number of measurements to gather to estimate center of transit time
    tc = 0
    told = zeros(tcmax) #holds timestamps 
    bold = zeros(tcmax) #holds impact parameter estimate

    xpos=0.0 ; ypos=0.0 ; zpos=0.0
    opos=0.0

    #Loop over each planet
    for j ∈ 2:nbody

        for i ∈ 1:npt
            time = sol.t[i] #time
            p[1:nvd2bd] = sol.u[i][1:nvd2bd] #momentum
            q[1:nvd2bd] = sol.u[i][nvd2bd + 1:nbody * NVEC] #position

            SR[1:nvecd2] = q[1:nvecd2] #store star position 
            SV[1:nvecd2] = p[1:nvecd2] ./ mass[1] #store star velocity

            # position of the planet relative to the star.
            xpos = (q[nvecd2 * j - 2] - SR[1])
            ypos = (q[nvecd2 * j - 1] - SR[2])
            zpos = (q[nvecd2 * j    ] - SR[3])

            # Rvec[1:nvecd2]=Array([x1, y1, z1]) #store planet position

            if xpos > 0.0 && opos <= 0.0
                tc = 0 #We are on the transit side.
            end
            if xpos <= 0.0
                tc = -1 #We are on the eclipse side, so skip and wait
            end

            if tc >= 0 #are on the transit side, so now we check for a transit
                b = zpos*zpos + ypos*ypos
                b = sqrt(b) / rstar
                if b < 1.1 #This really should be: b < 1+rprs.  Update when planet radius is added.  
                    tc += 1
                    told[tc] = time #record current time stamp
                    bold[tc] = b #record current position.
                end
            end

            if tc == tcmax #we need 5 transit points to do the calculation
                rt1 = (bold[4]-bold[5]) / (bold[2] - bold[3]) #compare slopes
                rt2 = (bold[3]-bold[4]) / (bold[1] - bold[2]) #compare slopes
                
                if rt1 < 0 || rt2 < 0 #see if we passed the minimum
                    tc = -1 #no point in searching for the minimum
                    s1 = (bold[1] - bold[2]) / (told[1] - told[2]) #Calculate slope
                    s2 = (bold[4] - bold[5]) / (told[4] - told[5])

                    b1 = bold[1] - s1*told[1]
                    b2 = bold[5] - s2*told[5]
                    if nTT[j] <= npt #prevent overflow
                        nTT[j] += 1
                        TT[j,nTT[j]] = (b2 -b1) / (s1 - s2) #estimate of mid transit time
                    end

                else

                    tc = tcmax-1 #no minimum, so delete earliest point

                    for k ∈ 1:tcmax-1
                        told[k] = told[k+1]
                        bold[k] = bold[k+1]
                    end
                end


            end

            #update old y-position
            opos = xpos
        end

        

    end

    nTT_max = maximum(nTT)
    TT = [TT[i, j] for i in 1:nbody, j in 1:nTT_max];

    return nTT, TT

end


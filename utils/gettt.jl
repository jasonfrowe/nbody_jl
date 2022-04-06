function getTT(sol, mass::AbstractVector, G::AbstractFloat)
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
    Rvec = zeros(nvecd2) #current position vector
    Rvec_old = zeros(nbody,nvecd2) #store reprevious position vector
    t_old = 0.0 #store previous time stamp

    for i ∈ 1:npt
        time = sol.t[i] #time
        p[1:nvd2bd] = sol.u[i][1:nvd2bd] #momentum
        q[1:nvd2bd] = sol.u[i][nvd2bd + 1:nbody * NVEC] #position

        SR[1:nvecd2] = q[1:nvecd2] #store star position 
        SV[1:nvecd2] = p[1:nvecd2] ./ mass[1] #store star velocity

        #now we get the position for each planet, skipping the star
        for j in 2:nbody
            x1 = (q[nvecd2 * j - 2] - SR[1])
            y1 = (q[nvecd2 * j - 1] - SR[2])
            z1 = (q[nvecd2 * j    ] - SR[3])

            Rvec[1:nvecd2]=Array([x1, y1, z1]) #Position of planet relative to star

            #check for transit
            if x1 > 0 && Rvec[2] / (Rvec_old[j, 2] + eps()) < 0 && i > 1
                
                ttrans = (t_old * Rvec[2] - time * Rvec_old[j,2]) / (Rvec[2] - Rvec_old[j,2])
                
                nTT[j] += 1
                TT[j, nTT[j]] = ttrans
                
                
                #Now we compute orbital solution

#                 RR=sum((Rvec).^2)
#                 R=sqrt(RR)
#                 Vvec=p[nvecd2*j-2:nvecd2*j]./mass[j] .- SV
#                 VV=sum((Vvec).^2)
#                 hh = ( (Rvec[2]*Vvec[3] - Rvec[3]*Vvec[2])^2 + (Rvec[3]*Vvec[1] - Rvec[1]*Vvec[3])^2 +
#                     (Rvec[1]*Vvec[2] - Rvec[2]*Vvec[1])^2 )
#                 gmpm=G*(mass[j]+mass[1])

#                 asemi1=(2/R - VV/gmpm)^(-1)
#                 fac=(1-hh/(gmpm*asemi1))
#                 if (fac)>0
#                     eccn1=sqrt(fac)
#                 else
#                     eccn1=eps()
#                 end

#                 th=asin(Rvec[2]/R)

#                 if eccn1<1.0 #only complete for bounded orbits
#                     cosf1 = 1.0/eccn1 * ( asemi1*(1.0-eccn1*eccn1)/R - 1.0 )
                    
#                     if abs(cosf1)<1
                        
#                         f=acos(cosf1)
#                         w=th-f
#                         cosf2 = cos(-w)

#                         #Now get the Eccentric Anomaly
#                         cosE1 = (cosf1+eccn1)/(1+eccn1*cosf1)
#                         cosE2 = (cosf2+eccn1)/(1+eccn1*cosf2)
#                         E1=acos(cosE1)
#                         E2=acos(cosE2)

#                         #and the orbital period
#                         Pd2pi=sqrt(asemi1^3/(G*(mass[j]+mass[1])))

#                         #Use the eccentric anomaly to predict when the mid-transit should of been
#                         t1=Pd2pi*(E1-eccn1*sin(E1))
#                         t2=Pd2pi*(E2-eccn1*sin(E2))
#                         dt=t2-t1 #This is the time correction
#                         ttrans=time+dt #this is the estimate of the centre of transit time

#                         nTT[j]+=1
#                         TT[j,nTT[j]]=ttrans
#                     end
#                 end

            end
            Rvec_old[j, 1:nvecd2] = Rvec[1:nvecd2]

        end
        t_old = time


    end

    nTT_max = maximum(nTT)
    TT = [TT[i, j] for i in 1:nbody, j in 1:nTT_max];
    
    return nTT, TT
    
end;
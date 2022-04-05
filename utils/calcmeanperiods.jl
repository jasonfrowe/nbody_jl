function calcmeanperiods(nTT,TT)

    #Requires pmodel

    nbody = size(nTT)[1]
    tt_periods = zeros(nbody)
    tt_T0 = zeros(nbody)
    for nb ∈ 2:nbody
        tt = TT[nb, 1:nTT[nb]]
        ntt = size(tt)[1]
        xp = zeros(2) #T0, Period
        xp[1] = minimum(tt) #guess for T0
        
        per = 0.0
        i = 0
        for j in 2:ntt
            per += tt[j] - tt[j-1]
            i += 1
        end
        per = per / i
        xp[2] = per #guess for periods

        tz = zeros(ntt)
        pfit = curve_fit(pmodel, tt, tz, xp; x_tol=1.e-12, g_tol=1.0e-14, maxIter=10000)
        tt_T0[nb] = pfit.param[1]
        tt_periods[nb] = pfit.param[2]
    end
    return tt_periods, tt_T0
end;

function pmodel(tt, xp)
    #Simple linear model used to estimate T0, period from transit times
   
    t0 = xp[1]
    per = xp[2]
    
    ph=(tt .- t0) ./ per
    omc=ph .- floor.(ph .+ 0.5)
    
    return omc
    
end;
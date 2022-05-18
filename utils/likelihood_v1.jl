function likelihood_v1(modelpars, extrapars)
    
    ep, tspan, G, nbody = extrapars

    mass = modelpars[1 : nbody]
    periods = modelpars[nbody + 1 : 2 * nbody]
    T0 = modelpars[2 * nbody + 1 : 3 * nbody]
    sqecosω = modelpars[3 * nbody + 1 : 4 * nbody]
    sqesinω = modelpars[4 * nbody + 1 : 5 * nbody]

    sol = calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)
    nTT, TT = getTT(sol, mass, G)
    #println(TT)

    large = 1.0e30 # a large number
    chi = 0.0
    for i ∈ 2:length(nTT_obs)

        if nTT[i] > 0
            for j ∈ 1:nTT_obs[i]

                chi += minimum( (((TT_obs[i,j] * DAYS) .- TT[i,1:nTT[i]]) ./ (TTerr_obs[i, j] * DAYS)).^2 )
                
            end
        else
            chi += large
        end
    end
    ll=-0.5*chi

    return ll

end;
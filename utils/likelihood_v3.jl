function likelihood(modelpars, extrapars)
    
    ep, tspan, G, nbody = extrapars

    mass = modelpars[1 : nbody]
    periods = modelpars[nbody + 1 : 2 * nbody]
    T0 = modelpars[2 * nbody + 1 : 3 * nbody]
    sqecosω = modelpars[3 * nbody + 1 : 4 * nbody]
    sqesinω = modelpars[4 * nbody + 1 : 5 * nbody]

    sol = calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)
    nTT, TT = getTT(sol, mass, G)
    #println(TT)

    ##Remove any slope from the TTVs -- so periods should be the TTV period from the data

    large = 1.0e30 # a large number
    chi = 0.0

    for i ∈ 2:length(nTT_obs)
        trnum_obs = Int.(floor.( (TT_obs[i,1:nTT_obs[i]] .* DAYS .- T0[i]) ./ periods[i] .+ 0.5 ))
        trnum = Int.( floor.( (TT[i,1:nTT[i]] .- T0[i]) ./ periods[i] .+ 0.5 ) )

        j=0
        for trnum1 ∈ trnum_obs
            j+=1
            pa = indexin(trnum1,trnum)
            if pa .== nothing
                chi += large
            else
                ## Update this to compare the OMC value at each numbered transit 
                chi += ((TT[i,pa[1]] - TT_obs[i,j] * DAYS) / (TTerr_obs[i,j] * DAYS))^2
            end

        end

        # println(trnum_obs)
        # println(trnum)

    end

    ll=-0.5*chi

    return ll

end;
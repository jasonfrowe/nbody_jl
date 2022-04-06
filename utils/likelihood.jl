function likelihood(mass, periods, T0, ep, sqecosω, sqesinω, tspan, G, nbody)

    sol = calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)
    nTT, TT = getTT(sol, mass, G)

    chi = 0.0
    for i ∈ 2:length(nTT_obs)

        for j ∈ 1:nTT_obs[i]

            chi += minimum( (((TT_obs[i,j] * DAYS) .- TT[i,1:nTT[i]]) ./ (TTerr_obs[i, j] * DAYS)).^2 )
            
        end
    end
    ll=-0.5*chi

    return ll

end;
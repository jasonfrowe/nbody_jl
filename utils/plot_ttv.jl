#function plot_TTV(chain,TTdf,periods,T0,nb)
function plot_TTV(chain, modelpars, extrapars, TT_obs, TTerr_obs, nTT_obs, nb)

    #Get parameters from input arguements 
    ep, tspan, G, nbody = extrapars

    periods_model = modelpars[nbody + 1 : 2 * nbody]
    T0_model = modelpars[2 * nbody + 1 : 3 * nbody]

    trace = []
    
    #nc = 5
    nc1 = Integer(floor(size(chain)[1]/2)) #this is equiv to burnin
    nchain = size(chain)[1]
    nsample = Integer(floor((nchain-nc1)/100))+1
    
    # We loop through the MCMC models and plot those.  
    i = 0
    for nc ∈ nc1:nsample:nchain
        
        i += 1

        mass = chain[nc][1][1 : nbody]
        periods = chain[nc][1][nbody + 1 : 2 * nbody]
        T0 = chain[nc][1][2 * nbody + 1 : 3 * nbody]
        sqecosω = chain[nc][1][3 * nbody + 1 : 4 * nbody]
        sqesinω = chain[nc][1][4 * nbody + 1 : 5 * nbody]
    
        sol = calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)
        # sol = calc_nbody(tstart, tend, chain[nc][1][1], chain[nc][1][2], 
        #     chain[nc][1][3], chain[nc][1][4], ep, dt)

        nTT, TT = getTT(sol, mass, G) #Get model transit times
        # print(nTT)

        P1 = periods_model[nb] / DAYS
        T1 = T0_model[nb]
        
        y = (TT[nb, 1:nTT[nb]] .- T1) ./ DAYS ./ (P1)
        y = (y .- floor.(y .+ 0.5)) .* (P1 * 24 * 60)
        # y=(y .- stat.median(y))
        x = TT[nb, 1:nTT[nb]] / DAYS
                
        # println(x)
        # println(y)
        if i == 1
            trace = scatter(x=x, y=y, mode="lines+markers", opacity=0.2)
        else
            trace1 = scatter(x=x, y=y, mode="lines+markers", opacity=0.2)
            trace = vcat(trace, trace1)
        end
    end
    
    y = (TT_obs[nb, 1:nTT_obs[nb]] .- T0_model[nb] / DAYS) ./ (periods_model[nb] / DAYS)
    y = (y .- floor.(y .+ 0.5)) .* periods_model[nb] ./DAYS *24 * 60
    
    x = TT_obs[nb, 1:nTT_obs[nb]]
    err = TTerr_obs[nb, 1:nTT_obs[nb]] .* (24.0 * 60.0)

    trace1 = scatter(x=x, y=y, mode="markers",
        error_y=attr(type="data", array=err, visible=true),marker=attr(color="black"))
    trace = vcat(trace, trace1)
    
    
    layout = Layout(;xaxis = attr(title="Time (days)"),
                     yaxis = attr(title="O-C (minutes)"),
                     showlegend = false,)
    
    
    plot(trace, layout)
    
    
end;
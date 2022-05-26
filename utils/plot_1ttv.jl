function plot_1ttv(nTT_obs, TT_obs, TTerr_obs, nTT, TT, modelpars, extrapars, nb)

    #Get parameters from input arguements 
    ep, tspan, G, nbody = extrapars

    periods_model = modelpars[nbody + 1 : 2 * nbody]
    T0_model = modelpars[2 * nbody + 1 : 3 * nbody]

    P1 = periods_model[nb] / DAYS
    T1 = T0_model[nb]
    y = (TT[nb, 1:nTT[nb]] .- T1) ./ DAYS ./ (P1)
    y = (y .- floor.(y .+ 0.5)) .* (P1 * 24 * 60)
    x = TT[nb, 1:nTT[nb]] / DAYS
    trace = scatter(x=x, y=y, mode="lines+markers", opacity=0.2)
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

end

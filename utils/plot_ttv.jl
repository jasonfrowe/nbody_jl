#function plot_TTV(chain,TTdf,periods,T0,nb)
function plot_TTV(chain, periods, T0, nTT_obs, TT_obs, TTerr_obs)

    trace = []
    
    #nc = 5
    nc1 = Integer(floor(size(chain)[1]/2)) #this is equiv to burnin
    nchain = size(chain)[1]
    nsample = Integer(floor((nchain-nc1)/100))+1
    
    i=0
    for nc in nc1:nsample:nchain
        
        i+=1
    
        sol=calc_nbody(tstart,tend,chain[nc][1][1],chain[nc][1][2],chain[nc][1][3],chain[nc][1][4],ep,dt)
        nTT,TT=getTT(sol, chain[nc][1][1], G)
        
        P1=periods[nb]/days
        T1=T0[nb]
        
        y=(TT[nb,1:nTT[nb]] .- T1) ./ days ./ (P1)
        y=(y.-floor.(y .+ 0.5)) .* (P1 *24 * 60)
        # y=(y .- stat.median(y))
        x=TT[nb,1:nTT[nb]]/days
                
        #println(x)
        #println(y)
        if i==1
            trace=scatter(x=x,y=y,mode="lines+markers",opacity=0.2)
        else
            trace1=scatter(x=x,y=y,mode="lines+markers",opacity=0.2)
            trace=vcat(trace,trace1)
        end
    end
    
    y=(TTdf[TTdf.planet_ind .== nb, :].tt .- T0[nb]/days) ./ (periods[nb]/days)
    y=(y.-floor.(y .+ 0.5)) .* periods[nb] ./days *24 * 60
    #y=(y .- stat.median(y))
    x=TTdf[TTdf.planet_ind .== nb, :].tt
    err=TTdf[TTdf.planet_ind .== nb, :].tt_err.*(24.0*60.0)
    
    trace1=scatter(x=x,y=y,mode="markers",
        error_y=attr(type="data", array=err, visible=true))
    trace=vcat(trace,trace1)
    
    
    layout = Layout(;xaxis=attr(title="Time (days)"),
                     yaxis=attr(title="O-C (minutes)"),
                     showlegend=false,)
    
    
    plot(trace, layout)
    
    
end;
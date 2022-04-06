function plotTTVs(tt_T0, tt_periods, nTT, TT)

    trace = []
    nbody = size(tt_T0)[1]
    for nb ∈ 2:nbody

        P1 = tt_periods[nb]
        T1 = tt_T0[nb]

        y = (TT[nb,1:nTT[nb]] .- T1) ./P1 
        y = (y .- floor.(y .+ 0.5)) * P1 ./ DAYS *24 * 60
        x = TT[nb, 1:nTT[nb]] / DAYS

        if nb == 2
            trace = scatter(x=x,y=y,mode="lines+markers",opacity=1)
        else
            trace1 = scatter(x=x,y=y,mode="lines+markers",opacity=1)
            trace = vcat(trace,trace1)
        end

    end

    layout = Layout(;xaxis=attr(title="Time (days)"),
                         yaxis=attr(title="O-C (minutes)"),
                         showlegend=false,width=800,height=600)


    plot(trace, layout)
end;

function plot_orbit(df,nbody,planet_names)
    
    ncolour=size(colors.phase)
    
    trace=[]
    plot_colours=[]
    npt=size(df.time)[1]
    nsample=Integer(floor(npt/1000)+1)
    
    amax=0
    for i in 1:nbody
        
        nc=Integer(floor((i-1)/nbody*ncolour+1))
        
        xpos_label=join(["xpos",string(i)])
        ypos_label=join(["ypos",string(i)])
        if i==1
            plot_colours=colors.phase[nc]
            trace = scatter(x=df[:,xpos_label][1:nsample:npt], y=df[:,ypos_label][1:nsample:npt],
                        mode="markers",marker=attr(color=colors.phase[nc],size=3))
        else
            plot_colours=vcat(plot_colours,colors.phase[nc])
            trace1 = scatter(x=df[:,xpos_label][1:nsample:npt], y=df[:,ypos_label][1:nsample:npt],
                            mode="markers",marker=attr(color=colors.phase[nc],size=3))
            trace=vcat(trace,trace1)
        end
        
        amax_x=maximum(df[:,xpos_label][1:nsample:npt])
        amax_y=maximum(df[:,ypos_label][1:nsample:npt])
        amax=maximum([amax,amax_x,amax_y])
        
    end
    
    x=[amax*0.9]
    y=[amax*0.9]
    for i in 2:nbody
        push!(x,amax*0.9)
        push!(y,amax*(0.9-(i-1)/nbody*0.35))
    end
    
    trace1 = scatter(x=x, y=y, text=planet_names, mode="markers+text", 
        textposition="right",marker=attr(color=plot_colours)) #size_max=60,
    trace=vcat(trace,trace1)
    
    
    layout = Layout(;xaxis=attr(title="X (au)"),
                     yaxis=attr(title="Y (au)"))
    plot(trace, layout)
    
end;

function plot_asemi(df,nbody)
    
    ncolour=size(colors.phase)
    
    trace=[]
    npt=size(df[:,"time"])[1]
    nsample=Integer(floor(npt/1000)+1)
    for i in 2:nbody
        
        nc=Integer(floor((i-1)/nbody*ncolour+1))
        
        asemi_label=join(["asemi",string(i)])
        
        if i==2
            trace= scatter(x=df[:,"time"][1:nsample:npt], y=df[:,asemi_label][1:nsample:npt],
                    mode="lines",marker=attr(color=colors.phase[nc]))
        else
            trace1=scatter(x=df[:,"time"][1:nsample:npt], y=df[:,asemi_label][1:nsample:npt],
                    mode="lines",marker=attr(color=colors.phase[nc]))
            trace=vcat(trace,trace1)
        end
                    
    end
    
    layout = Layout(;xaxis=attr(title="Time (years)"),
                     yaxis=attr(title="Semi-Major Axis (au)"))
    plot(trace, layout)
    
end;

function plot_eccn(df,nbody)
    
    ncolour=size(colors.phase)
    
    trace=[]
    npt=size(df[:,"time"])[1]
    nsample=Integer(floor(npt/1000)+1)
    for i in 2:nbody
        
        nc=Integer(floor((i-1)/nbody*ncolour+1))
        eccn_label=join(["eccn",string(i)])
        
        if i==2
            trace= scatter(x=df[:,"time"][1:nsample:npt], y=df[:,eccn_label][1:nsample:npt],
                    mode="lines",marker=attr(color=colors.phase[nc]))
        else
            trace1=scatter(x=df[:,"time"][1:nsample:npt], y=df[:,eccn_label][1:nsample:npt],
                    mode="lines",marker=attr(color=colors.phase[nc]))
            trace=vcat(trace,trace1)
        end
                    
    end
    
    layout = Layout(;xaxis=attr(title="Time (years)"),
                     yaxis=attr(title="Eccentricity"))
    plot(trace, layout)
    
end;

function plot_energy(df,nbody)
    
    ncolour=size(colors.phase)

    npt=size(df.time)[1]
    nsample=Integer(floor(npt/1000)+1)
    trace = scatter(x=df[:,"time"][1:nsample:npt], 
                    y=df[:,"TotalEnergy"][1:nsample:npt]/stat.median(df[:,"TotalEnergy"][1:nsample:npt]) .- 1,
                    mode="lines",marker=attr(color=colors.phase[1],size=3))
    
    layout = Layout(;xaxis=attr(title="Time (years)"),
                     yaxis=attr(title="Total Energy"))
    plot(trace, layout)
    
end;

function makeplots_v3(df,nbody,planet_names)
    p1 = plot_orbit(df,nbody,planet_names)
    p2 = plot_energy(df,nbody)
    p3 = plot_asemi(df,nbody)
    p4 = plot_eccn(df,nbody)
    p = [p1 p2; p3 ; p4]
    p.plot.layout["showlegend"] = false
    p.plot.layout["width"] = 800
    p.plot.layout["height"] = 1200
    p.plot.layout["vertical_spacing"] = 0.02 #This does nothing
    p
end;
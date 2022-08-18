function mhgmcmc(x, ex, llx, β, priors, nbuffer = 0, buffer = [0])
    
    corβ = 0.05 #controls vector jump amplitude
    rsamp = rand() #draw a random number to decide which sampler to use 
    
    x1 = x.* 1.0 #pre-allocate
    
    if nbuffer == 0 || rsamp < 0.5 #Metropolis-Hastings
    
        npars = size(x)[1]

        loop = 1
        while loop == 1
        
            n = Integer(floor(rand() * npars) + 1) #pick random number

            if priors[n, 5] > 0

                x1 = x * 1.0
                x1[n] += β[n] * randn()

                loop = 0 
            end

        end
        
    else #DE-MCMC
        
        n = -1
        n1 = Integer(floor(rand() * nbuffer) + 1) #pick 2 random numbers
        n2 = Integer(floor(rand() * nbuffer) + 1)
        
        xd1 = buffer[n1][1]
        xd2 = buffer[n2][1]
        
        vectorjump=xd1 .- xd2
        x1= x .+ (vectorjump .* corβ)
        
    end
        
    #pre-allocate
    llx1 = 0.0
    ac = [1,1]
    alpha = 0.0
    
    lp = logprior(x1, priors) #calculate prior to make sure model is valid
    #println("MCMC: $(n), $(lp), $(x1)")
    if lp > -Inf
        llx1 = likelihood(x1, ex) + lp
        alpha = minimum([exp(llx1-llx),1.0])
    else
        alpha = -1.0
    end
    
    u = rand()
    
    if u <= alpha #accept trial
        ac = [0, n] 
    else         #reject trial
        x1 = x .* 1.0
        llx1 = llx * 1.0
        ac = [1, n]
    end
    
    return llx1, x1, ac
    
end;
 
function genchain(x0, ex, β, priors, nsample, nbuffer = 0, buffer = [0])
   
    p = Progress(nsample; showspeed=true)

    llx=likelihood(x0, ex)
    chain=[[x0,[1,1]]] #first chain value
    llxp1 = 0.0 #pre-allocate for scope
    xp2=x0.*1.0
    ac = [1,1]
    for i in 1:nsample
        if nbuffer == 0
            llxp1,xp2,ac=mhgmcmc(x0, ex, llx, β, priors)
        else
            llxp1,xp2,ac=mhgmcmc(x0, ex, llx, β, priors, nbuffer, buffer)
        end
        llx = llxp1 * 1.0
        #println(ac)
        push!(chain,[xp2,ac])
        x0=xp2.*1.0
        ProgressMeter.next!(p; showvalues = [(:i,i), (:llxp1,llxp1)])
    end
    
    return chain
    
end;
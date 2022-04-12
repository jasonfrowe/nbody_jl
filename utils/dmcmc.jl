function mhgmcmc(x, ex, llx, β, priors ; nbuffer = 0, buffer = [0])
    
    corβ = 0.05 #controls vector jump amplitude
    rsamp = rand() #draw a random number to decide which sampler to use 
    
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
        
        xd1 = pars2x(buffer[n1][1], fitpars)
        xd2 = pars2x(buffer[n2][1], fitpars)
        
        vectorjump=xd1 .- xd2
        x1= x .+ (vectorjump .* corβ)
        
    end
        
    #pre-allocate
    llx1 = 0.0
    ac = [1,1]
    
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
        ac=[0, n] 
    else         #reject trial
        x1 = x .* 1.0
        llx1 = llx * 1.0
        ac = [1, n]
    end
    
    return llx1, x1, ac
    
end;
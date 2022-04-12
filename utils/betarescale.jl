function betarescale(x, ex, β, nsample, burnin; imax=20)
    
    alow = 0.22 #alow, ahigh define the acceptance rate range we want
    ahigh = 0.28
    
    delta = 0.01  #parameter controling how fast corscale changes - from Gregory 2011.
    
    npars = length(x)   #Number of parameters
    acorsub = zeros(npars) 
    nacor = zeros(npars)       #total number of accepted proposals
    nacorsub = zeros(npars)    #total number of accepted proposals immediately prior to rescaling
    npropp = zeros(npars)      #total number of proposals
    nproppsub = zeros(npars)   #total number of proposals immediately prior to rescaling
    acrate = zeros(npars)      #current acrate 
    corscale = ones(npars)
    
    #inital run
    chain = genchain(x, ex, β, priors, nsample) 
    nchain = size(chain)[1]
    
    #calcalate initial values of npropp and nacor 
    for i ∈ burnin:nchain
        j = Int(chain[i][2][2])           #get accept flag value
        npropp[j] += 1            #update total number of proposals
        nacor[j] += 1 - minimum([1, chain[i][2][1]]) #update total number of accepted proposals
    end
    
    #update x
    xin = chain[nchain][1]
    
    acrate = nacor ./ npropp #inital acceptance rate
    
    afix = ones(npars)  #afix is an integer flag to indicate which beta entries need to be updated
    for i ∈ 1:npars
        if acrate[i] < ahigh && acrate[i] > alow   #we strive for an acceptance rate between alow,ahigh
            afix[i] = 0    #afix=1 : update beta, afix=0 : do not update beta
        end
    end
    
    #We will iterate a maximum of imax times - avoid infinite loops
    icount = 0   #counter to track iterations
    while sum(afix) > 0
        icount += 1  #track number of iterations
        
        if icount > 1
            npropp = nproppsub .* 1
            nacor = nacorsub .* 1
        end
        nacorsub = nacorsub .* 0  #reset nacorsub counts for each loop
        nproppsub = nproppsub .* 0 #reset nproppsub counts for each loop
        
        #Make another chain starting with xin
        betain = β .* corscale   #New beta for Gibbs sampling 
        chain = genchain(xin, ex, β, priors, nsample)
        
        xin = chain[nchain][1] #Store current parameter state 
        
        for i ∈ burnin:nchain #scan through Markov-Chains and count number of states and acceptances 
            j = Int(chain[i][2][2]) 
            #if acrate[j]>ahigh or acrate[j]<alow: 
            npropp[j] += 1            #update total number of proposals
            nacor[j] += 1 - minimum([1, chain[i][2][1]]) #update total number of accepted proposals
            nproppsub[j] += 1            #Update current number of proposals
            nacorsub[j] += 1 - minimum([1, chain[i][2][1]]) #Update current number of accepted proposals
        end
            
        for i ∈ 1:npars  #calculate acceptance rates for each parameter that is to updated 
            #calculate current acrates
            if nproppsub[i] > 0 && npropp[i] - nproppsub[i] > 0
                acrate[i] = nacorsub[i] / nproppsub[i]

                #calculate acorsub
                acorsub[i] = (nacor[i] - nacorsub[i]) / (npropp[i] - nproppsub[i])

                if afix[i] > 0
                    #calculate corscale
                    corscale[i] = ( abs(corscale[i]
                        * ((acorsub[i]+delta)*0.75/(0.25*(1.0-acorsub[i]+delta)))^0.25 ) )
                end
            end
        end
        
        println("Current Acceptance: ",acrate) #report acceptance rates
        for i ∈ 1:npars  #check which parameters have achieved required acceptance rate
            if acrate[i] < ahigh && acrate[i] > alow
                afix[i] = 0
            end
        end
        
        if icount > imax   #if too many iterations, then we give up and exit
            afix = zeros(npars)
            println("Too many iterations: icount > imax")
        end
    end
    
    return corscale

end;
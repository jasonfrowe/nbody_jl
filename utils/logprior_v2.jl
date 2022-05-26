function logprior(modelpars, priors)

    lp = 1.0 #initialize log-prior 

    npars = size(priors)[1]
    for i ∈ 1:npars

        if priors[i, 5] == 1

            # Uniform distribution
            if modelpars[i] > priors[i,3] || modelpars[i] < priors[i,2]
                lp = -Inf
            else
                lp += 0.0
            end

        elseif priors[i, 5] == 2
           
            # Guassian prior
            lp +=  -( (modelpars[i] - priors[i,1]) / priors[1,4] )^2

        elseif priors[i, 5] == 3
           
            # Truncated prior
            if modelpars[i] > priors[i,3] || modelpars[i] < priors[i,2]
                lp = -Inf
            else
                lp += -( (modelpars[i] - priors[i,1]) / priors[1,4] )^2
            end

        end
        #println("lp: $(i) $(lp) $(priors[i, :])")
    end

    return lp

end
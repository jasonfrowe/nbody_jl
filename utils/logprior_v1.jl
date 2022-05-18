function logprior_v1(modelpars, priors)

    lp = 1.0 #initialize log-prior 

    npars = size(priors)[1]
    for i ∈ 1:npars

        if priors[i, 5] == 1

            # if modelpars[i] > priors[i,3] || modelpars[i] < priors[i,2]
            #     lp = 0.0
            # end
            lp += log( pdf(Uniform(priors[i,2], priors[i,3]), modelpars[i]) )

        elseif priors[i, 5] == 2

            lp += log( pdf(Normal(priors[i,1], priors[i,4]), modelpars[i]) )

        elseif priors[i, 5] == 3

            lp += log( pdf( truncated(Normal(priors[i,1], priors[i,4]), priors[i,2], priors[i,3]), modelpars[i]) )

        end
        #println("lp: $(i) $(lp) $(priors[i, :])")
    end

    return lp

end
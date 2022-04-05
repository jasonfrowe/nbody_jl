function readtt(nbody,ttdir,ttbase)

    nTT_obs = zeros(Int, nbody) #stores number of TTs for each lightcurve
    TT_obs = [] #stores the transit times
    TTerr_obs = [] #stores the measurement uncertainty in TTs
    push!(TT_obs, 0.0)
    push!(TTerr_obs, 0.0)
    for i ∈ 2:nbody
        #The filename with TTs - planet numbering starts at '1'
        filename=ttdir*ttbase*"."*lpad(i-1,2,"0")*".tt"
        #Check that the file exists
        if isfile(filename)

            #println("Found $(filename) $(isfile(filename)) ")
            
            #AbstractArray to hold the read in data
            time = [] #transit time 
            omc = [] #O-C value
            omc_err = [] #1 sigma error on O-C

            #Open the file
            f = open(filename)
            #go through the file line-by-line.  
            #This is our chance to do some data quality checks
            while eof(f) != true
                oneline = readline(f) #read in a single line
                onett = parse.(Float64, split(oneline)) #get the Vector 
                if onett[3] > 0.0 #the TTerr must be greater than zero.
                    push!(time, onett[1])
                    push!(omc, onett[2])
                    push!(omc_err, onett[3])
                end
            end
            close(f)
            #println(length(time))
            nTT_obs[i]=length(time)
            push!(TT_obs, time .+ omc)
            push!(TTerr_obs, omc_err)
        else
            println("Not found: $(filename) $(isfile(filename)) ")
        end
    end

    nTT_max=maximum(nTT_obs)
    TT_array=zeros(nbody,nTT_max)
    TTerr_array=zeros(nbody,nTT_max)
    for i ∈ 1:nbody
        for j ∈ 1:nTT_obs[i]
            TT_array[i,j]=TT_obs[i][j]
            TTerr_array[i,j]=TTerr_obs[i][j]
        end
    end

    return nTT_obs, TT_array, TTerr_array

end







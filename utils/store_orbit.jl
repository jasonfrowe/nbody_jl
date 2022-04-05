function store_orbit(sol::AbstractVector, 
    mass::AbstractVector, 
    nbody::AbstractVector, 
    NVEC::Integer)

    #Requires calxyzae.jl

    df = DataFrame()
    npt = size(sol.t)[1]
    nbodyvecd2 = Integer(nbody * NVEC / 2)
    nvecd2 = Integer(NVEC / 2)
    
    #Variables for threading
    # nthread=Threads.nthreads()
    
    df.time = sol.t
    
    for i ∈ 1:nbody
        
        x = zeros(npt)
        y = zeros(npt)
        z = zeros(npt)
        eccn = zeros(npt)
        asemi = zeros(npt)
    
        Threads.@threads for j ∈ 1:npt
            x[j], y[j], z[j], eccn[j], asemi[j] = calxyzae(i, sol.u[j][1:nbodyvecd2], 
               sol.u[j][nbodyvecd2+1:nbody*NVEC], mass, nbody, nvecd2)
        end

        df.t1 = x
        df.t2 = y
        df.t3 = z
        df.t4 = eccn
        df.t5 = asemi
        
        colname = join(["xpos", string(i)])
        rename!(df, :t1 => colname )

        colname = join(["ypos", string(i)])
        rename!(df, :t2 => colname )
    
        colname = join(["zpos", string(i)])
        rename!(df, :t3 => colname )
        
        colname = join(["eccn", string(i)])
        rename!(df, :t4 => colname )

        colname = join(["asemi", string(i)])
        rename!(df, :t5 => colname )
    
    end
    
    TotalEnergy = zeros(npt)
    t = 0
    Threads.@threads for i ∈ 1:npt
        TotalEnergy[i] = H_simd(sol.u[i][1:nbodyvecd2], sol.u[i][nbodyvecd2+1:nbody*NVEC], mass, t)
    end
    df.TotalEnergy = TotalEnergy
    
    return df
end
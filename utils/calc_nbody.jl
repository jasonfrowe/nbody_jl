function calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)

    #Convert from Keplerian to Cartesian state
    p0, q0 = get_ystart_v5(mass, periods, T0, ep, sqecosω, sqesinω);

    q0=SVector{size(q0)[1]}(q0); #store as Static Array for better speed
    p0=SVector{size(p0)[1]}(p0);
    mass=SVector{size(mass)[1]}(mass);
    param=mass;

    prob = HamiltonianProblem(H_simd, p0, q0, tspan, param); #set up Hamiltonian to solve

    #Calculate Model @jasonfrowe Should move abstol and relton to parameter file
    sol = solve(prob, KahanLi8(), abstol=1e-10, reltol=1e-10, dt=Δt);
    
    # #Store model
    # df = store_orbit(sol, mass, nbody, NVEC);
    #Extract TTs
    nTT,TT=getTT(sol, mass, G);

    #Calcalate Mean Periods
    try
        tt_period,tt_T0=calcmeanperiods(nTT,TT);
    catch
        tt_period = periods .* 1.0
        tt_T0 = T0 .* 1.0
    end

    #Get period correction
    periods_cor = periods .+ (periods .- tt_period)
    #T0_cor = T0 .+ (T0 .- tt_T0)
    T0_cor = TT[:,1]

    #redo model...
    p0, q0 = get_ystart_v5(mass, periods_cor, T0_cor, ep, sqecosω, sqesinω);

    q0=SVector{size(q0)[1]}(q0); #store as Static Array for better speed
    p0=SVector{size(p0)[1]}(p0);
    mass=SVector{size(mass)[1]}(mass);
    param=mass;

    prob = HamiltonianProblem(H_simd, p0, q0, tspan, param); #set up Hamiltonian to solve

    #Calculate Model @jasonfrowe Should move abstol and relton to parameter file
    sol = solve(prob, KahanLi8(), abstol=1e-10, reltol=1e-10, dt=Δt);   

    return sol

end
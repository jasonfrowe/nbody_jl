function calc_nbody(mass, periods, T0, ep, sqecosω, sqesinω, tspan)

    #Convert from Keplerian to Cartesian state
    p0, q0 = get_ystart_v4(mass, periods, T0, ep, sqecosω, sqesinω);

    q0=SVector{size(q0)[1]}(q0); #store as Static Array for better speed
    p0=SVector{size(p0)[1]}(p0);
    mass=SVector{size(mass)[1]}(mass);
    param=mass;

    prob = HamiltonianProblem(H_simd, p0, q0, tspan, param); #set up Hamiltonian to solve

    #Calculate Model @jasonfrowe Should move abstol and relton to parameter file
    sol = solve(prob, KahanLi8(), abstol=1e-10, reltol=1e-10, dt=Δt);

    return sol

end
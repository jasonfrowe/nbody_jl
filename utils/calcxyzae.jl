function calxyzae(i::Integer, 
    p::AbstractVector, 
    q::AbstractVector, 
    mass::AbstractVector, 
    nbody::AbstractVector, 
    nvecd2::Integer)

    CoM = zeros(Integer(nvecd2)) #center of mass
    CoV = zeros(Integer(nvecd2)) #center of mass
    for k in 1:nbody
        CoM += q[nvecd2 * k - 2:nvecd2 * k] * mass[k]
        CoV += p[nvecd2 * k - 2:nvecd2 * k]
    end
    CoM = CoM ./ sum(mass)
    CoV = CoV ./ sum(mass)
    x1 = (q[nvecd2 * i - 2] - CoM[1])
    y1 = (q[nvecd2 * i - 1] - CoM[2])
    z1 = (q[nvecd2 * i    ] - CoM[3])
    
    Rvec = Array([x1, y1, z1])
    RR = sum((Rvec).^2)
    R = sqrt(RR)
    Vvec = p[nvecd2 * i - 2:nvecd2 * i] ./ mass[i] .- CoV
    VV = sum((Vvec).^2)
    hh = ( (Rvec[2] * Vvec[3] - Rvec[3] * Vvec[2])^2 + (Rvec[3] * Vvec[1] - 
          Rvec[1] * Vvec[3])^2 + (Rvec[1] * Vvec[2] - Rvec[2] * Vvec[1])^2 )
    gmpm = G * (mass[i] + mass[1])

    asemi1 = (2 / R - VV / gmpm)^(-1)
    fac = (1 - hh / (gmpm * asemi1))
    if (fac) > 0
        eccn1 = sqrt(fac)
    else
        eccn1 = 0.0
    end
    
    # Convert to MKS
    x1 = x1 / AU
    y1 = y1 / AU
    z1 = z1 / AU
    asemi1 = asemi1 / AU

    return x1,y1,z1,eccn1,asemi1
end;

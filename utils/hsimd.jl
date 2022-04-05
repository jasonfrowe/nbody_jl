function H_simd(p, q, mass, t)
    
    nbody = length(mass)
    nvecd2 = Integer(NVEC / 2)
    
    h1 = 0.0
    h2 = 0.0
    for i ∈ 1:nbody
        pp = 0
        @simd for k ∈ 0:nvecd2 - 1
            pp += p[i * nvecd2 - k] * p[i * nvecd2 - k]
        end
        h1 += pp / (2 * mass[i])   
        
        for j ∈ i + 1:nbody
            qq = 0.0
            @simd for k ∈ 0:nvecd2-1
                qq += (q[i * nvecd2 - k] - q[j * nvecd2 - k]) * (q[i * nvecd2 - k] - q[j * nvecd2 - k])
            end
            h2 += mass[i] * mass[j] / sqrt(qq)
        end
    end
    
    h = h1 - G * h2
    return h
end;
"""
    function guide(x, P::SIRguided, h, z)

    evolve the guided process for one time step from 
    present state x using htransform h with randomness z
    
    returns new state and logweight 
"""
function guide!(xnext, xcurrent, P::SIRguided, h, z, infected)
    @assert length(xnext)==length(xcurrent)==length(h)==length(z)
    logweight = 0.0 
    for i ∈ eachindex(xcurrent)
        if xcurrent[i]==_S_
            ni = nr_infected_neighb(xcurrent, P.𝒩, i)
            p = pS(P.λ * ni * P.τ) .* h[i] # allocates 10

            ñi = infected[i]
            p̃ = pS(P.λ * ñi * P.τ) .* h[i]
            logweight += log(sum(p)) - log(sum(p̃))
        elseif xcurrent[i]==_I_
            p = pI(P.μ * P.τ) .* h[i]
        elseif xcurrent[i]==_R_
            p = pR(P.ν * P.τ) .* h[i]
        end
        xnext[i] = rand𝒳(z[i], p/sum(p)) 
    end
    logweight
end


"""
    function forward(P::SIRguided, Π, B, Z, logw)

    simulate guided process using prior Π on the initial state (indexed by "1")
 
    B contains output of backward filter (contains n_times vectors, where each of these vectors
    contains n_particle vectors in ℝ³)

    Z contains innovations (random numbers for simulating the guided process)

    returns simulated path and loglikelihood
"""
function forward(P::SIRguided, Π, B, Z, logw)
    n_steps, n_particles = length(Z), length(Π)

    # sample initial state
    X = Vector{State}(undef, n_particles)
    ll = logw
    
    for i in 1:n_particles
        p = Π[i] .* B[1][i]
        X[i] = rand𝒳(Z[1][i], p/sum(p))
        ll += log(sum(p))
    end
    

    Xs = [copy(X)]
    for t in 2:n_steps
        lw = guide!(X, Xs[t-1], P, B[t], Z[t], P.ℐ[t-1])
        ll += lw
        push!(Xs, copy(X))
    end
    Xs, ll
end


forward(P, Π, B, logw) = (Z) -> forward(P, Π, B, Z, logw)

# in place version
function forward!(Xs, P::SIRguided, Π, B, Z, logw)
    n_steps, n_particles = length(Z), length(Π)

    # sample initial state
    #X = Vector{State}(undef, n_particles)
    ll = logw
    
    for i in 1:n_particles
        p = Π[i] .* B[1][i]
        Xs[1][i] = rand𝒳(Z[1][i], p/sum(p))
        ll += log(sum(p))
    end
    


    for t in 2:n_steps
        lw = guide!(Xs[t], Xs[t-1], P, B[t], Z[t], P.ℐ[t-1])
        ll += lw
    end
    ll
end



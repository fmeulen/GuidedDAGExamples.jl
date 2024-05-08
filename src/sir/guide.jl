"""
    function guide(x, P::SIRguided, h, z)

    evolve the guided process for one time step from 
    present state x using htransform h with randomness z
    
    returns new state and logweight 
"""
function guide!(x, P::SIRguided, h, z, infected)
    @assert length(x)==length(h)==length(z)
    logweight = 0.0 
    for i ∈ eachindex(x)
        if x[i]==_S_
            ni = nr_infected_neighb(x, P.𝒩, i)
            p = pS(P.λ * ni * P.τ) .* h[i]

            ñi = infected[i]
            p̃ = pS(P.λ * ñi * P.τ) .* h[i]
            logweight += log(sum(p)) - log(sum(p̃))
        elseif x[i]==_I_
            p = pI(P.μ * P.τ) .* h[i]
        elseif x[i]==_R_
            p = pR(P.ν * P.τ) .* h[i]
        end
        x[i] = rand𝒳(z[i], p/sum(p)) 
    end
    logweight
end


"""
    function forward(P::SIRguided, Π, B, Z)

    simulate guided process using prior Π on the initial state (indexed by "1")
 
    B contains output of backward filter (contains n_times vectors, where each of these vectors
    contains n_particle vectors in ℝ³)

    Z contains innovations (random numbers for simulating the guided process)

    returns simulated path and loglikelihood
"""
function forward(P::SIRguided, Π, B, Z)
    n_steps, n_particles = length(Z), length(Π)

    # sample initial state
    X = Vector{State}(undef, n_particles)
    z = Z[1]
    ll = 0.0
    for i in 1:n_particles
        p = Π[i] .* B[1][i]
        X[i] = rand𝒳(z[i], p/sum(p))
        ll += log(sum(p))
    end
    

    Xs = [deepcopy(X)]
    for t in 2:n_steps
        lw = guide!(X, P, B[t], Z[t], P.ℐ[t-1])
        ll += lw
        push!(Xs, deepcopy(X))
    end
    Xs, ll
end

forward(P, Π, B) = (Z) -> forward(P, Π, B, Z)


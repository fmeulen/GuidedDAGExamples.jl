"""
    function guide(x, P::SIRguided, h, z, infected_auxiliary)

    evolve the guided process for one time step from 
    present state x using htransform h with randomness z
    
    infected_auxiliary: nr of infected under auxiliary process

    returns new state and logweight 
"""
function guide!(xnext, xcurrent, P::SIRguided, h, z, infected_auxiliary)
    @assert length(xnext)==length(xcurrent)==length(h)==length(z)
    logweight = 0.0 
    for i ∈ eachindex(xcurrent)
        if xcurrent[i]==_S_
            ni = nr_infected_neighb(xcurrent, P.𝒩, i)
            p = pS(P.λ * ni * P.τ) .* h[i] # allocates 10

            # following three lines should not be part of this function
            ñi = infected_auxiliary[i]
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
function forward(P::SIRguided, Π, B, Z, l0)
    n_steps, n_particles = length(Z), length(Π)

    ll = l0
    
    # sample initial state
    X = Vector{State}(undef, n_particles)
    for i in 1:n_particles
        p = Π[i] .* B[1][i]
        X[i] = rand𝒳(Z[1][i], p/sum(p))
        ll += log(sum(p)) # contribution to loglik of pullback of prior Π
    end

    Xs = [copy(X)]
    for t in 2:n_steps
        lw = guide!(X, Xs[t-1], P, B[t], Z[t], P.ℐ[t-1])
        ll += lw
        push!(Xs, copy(X))
    end
    # compute contribution of fully simulated guided path to likelihood, using B and (perhaps) 𝒪
    Xs, ll
end

# separate function to compute the logweight
"""
    logweight(Xs, Π, B, 𝒪, O)

    Xs: simulated guided process
    Π: prior on initial state
    B: backward filter
    𝒪: observations
    O: emission matrix to observations (assumed to be constant for all observations)
"""
function loglikelihood(Xs, Π, B, 𝒪, O)
# compute logweight for guided path Xs that was obtained using backward filter B
    ll = 0.0
    # pullback from prior
    for i in 1:n_particles
        p = Π[i] .* B[1][i]
        ll += log(sum(p)) # contribution to loglik of pullback of prior Π
    end
    # weights at all time steps exclusing contribution of observations
    n_times = length(𝒪)
    for t in 1:n_times
        xt = Xs[t]
        bt = B[t]
        
        for i in 1:n_particles
            x = xt[i] 
            g1 = bt[i]
            g2 = t==n_times ?  SA_F64[1, 1, 1] : B[t+1][i]
            if x ==_S_
                ni = nr_infected_neighb(xt, P.𝒩, i)
                p = pS(P.λ * ni * P.τ)  
            elseif x==_I_
                p = pI(P.μ * P.τ) 
            elseif x==_R_
                p = pR(P.ν * P.τ) 
            end
            ll += (t ≠ n_times)*log(dot(p,g2)) - log(g1[ind(x)]) 
        end    
    end
    # contribution from observations
    for t in eachindex(𝒪)
        for k in 𝒪[t].ind
            observation = 𝒪[t].x[k]
            xtk = Xs[t][k]
            ll += log(O[ ind(xtk) , ind(observation)]) # ind is a function that converts S,I,R to 1,2,3
        end
    end
    ll        
end



# in place version
function forward!(X, P::SIRguided, Π, B, Z, l0)
    n_steps, n_particles = length(Z), length(Π)
    ll = l0

    # sample initial state
    for i in 1:n_particles
        p = Π[i] .* B[1][i]
        X[1][i] = rand𝒳(Z[1][i], p/sum(p))
        ll += log(sum(p))
    end
    
    for t in 2:n_steps
        lw = guide!(X[t], X[t-1], P, B[t], Z[t], P.ℐ[t-1])     
        ll += lw
    end
    ll
end



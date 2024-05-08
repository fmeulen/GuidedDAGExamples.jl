######################## backward filtering #########################
# h is a vector of Svectors of length 3
# 
function fuse!(O::Observation, h)
    id = O.ind
    println(id)
    for i in eachindex(id)
        k = id[i]  # this is an index in h that needs to be updated because we observe it
        h[k] = O.h[i] .* h[k]
    end
end

# back kernel for one individual
κ̃(P::SIRguided ,ninfected::Number) = hcat(pS(P.λ * P.τ * ninfected), pI(P.μ*P.τ), pR(P.ν*P.τ))'

"""
    h hfun at time t+1
    n vector of infected individuals at time t
"""
function pullback!(h, ninfected, P::SIRguided) 
    for i in eachindex(h)    
        h[i] = κ̃(P, ninfected[i]) * h[i]
    end
end

function normalize!(h)    
    for i in eachindex(h)
        h[i] = h[i]/sum(h[i])
    end
end

"""
    nr_infected_neighb(x,𝒩,i)

    Computes number of infected neighbours for the i-th individual in configuration x (at a particular time)
    If x[i] !== _S_ then it is set to zero (because not needed)

    𝒩 = set_neighbours(8,2)
    X = [_S_, _I_, _S_, _I_, _R_, _S_, _I_, _I_]
    for i in 1:8
        @show  nr_infected_neighb(X, 𝒩, i)
    end

"""
nr_infected_neighb(x, 𝒩, i) = x[i] == _S_ ? sum(x[𝒩[i]].==_I_) : 0

"""
count_infections_at_t(x, 𝒩)

    count at one time instance for one particle
"""
count_infections_at_t(x, 𝒩) =[nr_infected_neighb(x, 𝒩, i) for i in eachindex(x)]  


""" count_infections(XX, 𝒩)

    n_particles = 8
    𝒩 = set_neighbours(8)
    X1 = [_S_, _I_, _S_, _I_, _R_, _S_, _I_, _I_]
    for i in 1:8
        @show  nr_infected_neighb(X1, 𝒩, i)
    end

    X2 = [_I_, _L_, _S_, _L_, _L_, _S_, _I_, _I_]
    XX = [X1, X2]
    count_infections(XX, 𝒩)
"""
count_infections(X, 𝒩) = [count_infections_at_t(x, 𝒩)  for x ∈ X]





function backward(P::SIRguided, 𝒪)
    n_times = length(𝒪)
    n_particles = length(𝒪[1].x)
    
    h = [SA_F64[1, 1, 1]  for _ in 1:n_particles]
    fuse!(𝒪[n_times], h)
    hs = [copy(h)]
    for t in n_times-1:-1:1
        pullback!(h, P.ℐ[t], P)
        
        fuse!(𝒪[t], h)
        normalize!(h)
        pushfirst!(hs, copy(h))
    end
    hs
end


######################## backward filtering #########################

"""
    fuse!(O::Observation, h)

    O: contains one observation
    h: is a vector of Svectors of length 3

    those indices in h where there is an observatoin in O are fused
"""
function fuse!(O::Observation, h)
    ids = O.ind
    for i in eachindex(ids)
        k = ids[i]  # this is an index in h that needs to be updated because we observe it
        h[k] = O.h[i] .* h[k]
    end
end

# backward kernel for one individual
κ̃(P::SIRguided ,ninfected::Number) = hcat(pS(P.λ * P.τ * ninfected), pI(P.μ*P.τ), pR(P.ν*P.τ))'

"""
    pullback!(h, ninfected, P::SIRguided)

    idea:
    h is hfun at time t+1
    ninfected is the  vector of infected individuals at time t
"""
function pullback!(h, ninfected, P::SIRguided)
    for i in eachindex(h)
        h[i] = κ̃(P, ninfected[i]) * h[i]
    end
end

"""
    normalize!(h)

    h: is a vector of Svectors of length 3
    each element gets mutated such that its elements sum to one
    sum of all log normalisation constants is returned
"""
function normalize!(h)
    s = 0.0
    for i in eachindex(h)
        si = sum(h[i])
        h[i] = h[i]/si
        s += log(si)
    end
    s
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

    (function does not allocate)
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


"""
    backward(P::SIRguided, 𝒪)

    Returns backward information filter using diagonalisation

    Additionally returns log likelihood contribution induced
    by normalisation (backw filtering to root node using prior Π is not
    included)
"""
function backward(P::SIRguided, 𝒪)
    n_times = length(𝒪)
    n_particles = length(𝒪[1].x)
    ll0 = 0.0

    h = fill(SA_F64[1, 1, 1], n_particles) #start with uniforms (should this be normalised?)
    fuse!(𝒪[n_times], h)                 # fuse in obs at time T ( == n_times)
    B = [copy(h)]                        # save guiding h in B array
    for t in n_times-1:-1:1              # for t = T-1 back to 1 do:
        pullback!(h, P.ℐ[t], P)          # pullback with h_{t+1} and 'known' nr. Infected at time t
        fuse!(𝒪[t], h)                   # fuse pullback h_t with obs_t

        lw = normalize!(h)               # normalise
        ll0 += lw                        # account for normalisation
        pushfirst!(B, copy(h))           # save guiding h in Barray
    end
    B, ll0
end



######################## backward filtering #########################
# h is a vector of Svectors of length 3
# 
function fuse(O::Observation, h)
    id = O.ind
    for i in eachindex(id)
        h[i] = O.h[i] .* h[id[i]]
    end
    h
end

# back kernel for one individual
κ̃(θ,ninfected,τ) = hcat(pS(θ.λ * τ * ninfected), pI(θ.μ*τ), pR(θ.ν*τ))'

"""
    h hfun at time i+1
    n vector of infected individuals at time i
"""
pullback(h, ninfected, θ, τ) = [κ̃(θ, ninfected[i], τ) * h[i] for i ∈ eachindex(h)]

# initalisation
# ninfected = [rand(Poisson(1),n_particles) for _ in 1:n_times]
# θ = (λ=1.0, ν =2.0, μ=0.5)
# τ = 0.1


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
count_infections_at_t(X, 𝒩) =[nr_infected_neighb(X, 𝒩, i) for i in eachindex(X)]

count_infections(XX, 𝒩) = float.([count_infections_at_t(x, 𝒩)  for x ∈ XX])


function backward(P::SIRguided, 𝒪, infected_neighbours)
    n_times = length(𝒪)
    n_particles = length(𝒪[1].x)
    θ = params(P)

    h_ = [SA_F64[1, 1, 1]  for _ in 1:n_particles]
    h = fuse(𝒪[n_times], h_)
    hs = [h]
    for t in n_times-1:-1:1
        h_ = pullback(h, infected_neighbours[t], θ, P.τ)
        h = fuse(𝒪[t], h_)
        pushfirst!(hs, copy(h))
    end
    hs
end

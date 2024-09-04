"""
    This function will have to:
        1. initialise all static CPD tables (in last year's code this was done in FactorisedMarkovChain)
        2. perform EP backward filtering given the observations
    I will resuse alot of code from earlier implementation
"""
function backwardEP(P::SIRguided, 𝒪)
    # need this info to make old code work without many edits
    N = n_particles
    T = n_times

    size_neighbourhood = length(P.𝒩[1])
# x0 is not passed to the function, I would not expect it here
    root = x0  


    # the following code is taken from SIR.jl
    # https://github.com/fmeulen/BackwardFilteringSIR/blob/main/SIR.jl

    δ = 0.0 # artificial par to become infected without infected neighbours
    τ = P.τ # discretisation time step of SIR model

    statespace = Dict([i => E for i in 1:N])
    parents = Dict(i => intersect(i-size_neighbourhood:i+size_neighbourhood, 1:N) for i in 1:N)
    
    # Parametric description of the entire forward model (I wish to pass δ and τ as well)
    SIR(θ, δ, τ) = FactorisedMarkovChain(statespace, parents, dynamics(θ, δ, τ), root, (N, T))

    # Instantiation with the true dynamics
    θ = [P.λ, P.μ, P.ν]
    G = SIR(θ, δ, τ)

    # Below is how I had implemented all required observations info
    obsparents = merge([Dict((i,t) => (i,t) for i in 𝒪[t].ind) for t in 1:T]...)
    obscpds = merge([Dict((i,t) => O for i in 𝒪[t].ind) for t in 1:T]...)
    obsstates = merge([Dict((i,t) => Int(𝒪[t].x[i]) for i in 𝒪[t].ind) for t in 1:T]...)

    obs = (obsparents, obscpds, obsstates)

# We don't want this contribution, can we get rid of it?    
    # Prior is required for contribution from root.
    Πroot =  Dict(i => [0.98, 0.02, 0.00] for i in 1:N)

    # Backward filter
    propagation = boyenkoller
    ms, logh =  backwardfiltering(G, propagation, false, obs, Πroot, size_neighbourhood)

    # now ms exists only for time steps 2:T, this is because in the past we had 'fixed root' x_1
    # we also want to compute the integral over prior x1 (what is called here logh) in a seperate function
    # so we just return only the guiding functions (what you call B) which we can get from messages ms
    # recall that we do not 'care' or concern ourselves with normalisation in the bw filter procedure

    # recall that in my implementation ms[t].htransforms == ms[t+1].approximatepullback, thus
    B = [[ms[2].approximatepullback[i] for i in 1:N], [[ms[t].factoredhtransform[i] for i in 1:N] for t in 2:T]...]

    B, logh
end

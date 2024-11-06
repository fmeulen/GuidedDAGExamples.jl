using Random

"""
    resample(particles, weights)

Resamples particles according to their weights using multinomial resampling.
Returns a new set of particles with duplicates based on their importance.

# Arguments:
- `particles`: A vector of particles.
- `weights`: A vector of weights, should sum to 1.

# Returns:
- `resampled_particles`: A new vector of resampled particles.


# Example particles and weights
particles = ["p1", "p2", "p3", "p4"]
weights = [0.1, 0.2, 0.4, 0.3]

# Perform resampling
new_particles = resample(particles, weights)

println("Resampled Particles: ", new_particles)

"""
function resample(particles, weights)
    N = length(particles)
    # Normalize the weights to sum to 1 (if not already normalized)
    weights /= sum(weights)
    
    # Generate N random indices based on the multinomial distribution
    resampled_indices = rand(Categorical(weights), N)
    
    # Resample particles based on the generated indices
    resampled_particles = particles[resampled_indices]
    
    return resampled_particles
end

"""
    resamp(lw, ess_threshold)

Function to exponentiate and normalize log-weights
lw = [log(1.2), log(0.3), log(10.1), log(0.4)]
resamp(lw, .2)
"""
function resamp(lw, ess_threshold)
    w = exp.(lw .- maximum(lw))
    println(w)
    nw = w/sum(w)
    ess = 1 / sum(nw.^2)
    println("ESS: $ess")
    #println(ess)
    if ess < ess_threshold
        println("ESS is below threshold, resampling takes place\n")
        return resample(1:length(lw), nw) #SMC.resample(nw)
    else
        return 1:length(lw)
    end
end


function move!((Z, X, ll), NR_MOVE_STEPS, P, Π, Zᵒ, Xᵒ, B, prior, blocks, δ, printskip)
    i = 1
    #ws = [ll]
    for i = 1:NR_MOVE_STEPS
        ll, Z, Zᵒ, X, Xᵒ, acc_ = mcmc_iteration!(P, Π, Z, Zᵒ, X, Xᵒ, B, ll, prior, blocks, δ, i, printskip)
    #    push!(ws, ll)
    end
     (Z=Z, X=X, logweight=ll)
end

function inititalise_particle(P, Π, B, prior)
    n_times, n_particles = length(B), length(B[1])
    Z = innovations(n_times, n_particles)
    X, ll  = forward(P, Π, B, Z, prior) 
    (Z=Z, X=X, logweight=ll)
end


"""
    smc(NR_SMC_STEPS, NUMPARTICLES, NR_MOVE_STEPS, P, Π, prior, blocks, δ, printskip)

NUMPARTICLES = 10 
NR_SMC_STEPS = 50
NR_MOVE_STEPS =50
printskip = 1000
δ = 0.01

out = smc(NR_SMC_STEPS, NUMPARTICLES, NR_MOVE_STEPS, P, Π, prior, blocks, δ, printskip)
"""
function smc(NR_SMC_STEPS, NUMPARTICLES, NR_MOVE_STEPS, P, Π, prior, blocks, δstep, printskip)
   # initialise particles 
    B = backward(P)
    particles = [inititalise_particle(P, Π, B, prior) for _ in 1:NUMPARTICLES]
    
    Zᵒ = [deepcopy(particles[1].Z) for _ in 1:NUMPARTICLES]
    Xᵒ = [deepcopy(particles[1].X) for _ in 1:NUMPARTICLES]

    ess_threshold = round(NUMPARTICLES/2; digits=0)    
    lls = []

    for j in 1:NR_SMC_STEPS
        println(j)
        # resample step
        log_weights =  getindex.(particles, :logweight)
        indices = resamp(log_weights, ess_threshold) 
        println(indices, "\n")
        particles = [particles[k] for k in indices]
        # move step (pCN)
        particles = map(i -> move!(particles[i], NR_MOVE_STEPS, P, Π, Zᵒ[i], Xᵒ[i], B, prior, blocks, δstep, printskip), 1:NUMPARTICLES)
        push!(lls, [particles[i].logweight for i in eachindex(particles)])
    end
    (particles=particles, logweights=lls)
end


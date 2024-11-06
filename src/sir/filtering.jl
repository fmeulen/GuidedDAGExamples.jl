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

function ess(lw) # lw = logweights
    w = exp.(lw .- maximum(lw))
    nw = w/sum(w)
    ess_ = 1.0 / sum(nw.^2)
    ess_, nw 
end
"""
    resamp(lw, ess_threshold)

Function to exponentiate and normalize log-weights
lw = [log(1.2), log(0.3), log(10.1), log(0.4)]
resamp(lw, .2)
"""
function resamp(lw, ess_threshold)
    ess_, nw = ess(lw)
    println("ESS: $ess_")
    if ess_ < ess_threshold
        println("ESS is below threshold, resampling takes place\n")
        return resample(1:length(lw), nw) #SMC.resample(nw)
    else
        return 1:length(lw)
    end
end


function construct_Ploc(P, i0, iT)
    𝒪loc = P.𝒪[i0:iT]
    ℐloc =initialise_infected_neighbours(𝒪loc) 
    SIRguided(P.ξ, P.λ, P.μ, P.ν, P.τ, P.𝒩, ℐloc, 𝒪loc, P.O) 
end

function convert_state_to_Π(x::State)
    if x == _S_
        return(@SVector([1,0,0]))
    elseif x == _I_
        return(@SVector([0,1,0]))
    else
        return(@SVector([0,0,1]))
    end
end

"""
guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble)

For all ensemble members (n_ensemble in total)
-   Sample guided process in {i0, i0+1, ... , iT}
-   Compute logweight of this path.
where Πloc is a prior on the state at time i0 for each ensemble member. 
prior: prior on parameters in P (μ, λ, ν)
"""
function guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble)
    Ploc = construct_Ploc(P,i0,iT)
    Bloc = backward(Ploc)

    n_times = ntimes(Ploc) #iT - i0 + 1
    n_particles = nparticles(Ploc) #length(P.𝒩)
    Z = innovations(n_times, n_particles)
    Xᵒ, logweight  = forward(Ploc, Πloc[1], Bloc, Z, prior)

    Xs = [Xᵒ]
    logweights = [logweight]
    
    # repeat all  other particles
    for i in 2:n_ensemble 
        Z = innovations(n_times, n_particles)
        Xᵒ, logweight  = forward(Ploc, Πloc[i], Bloc, Z, prior)
        push!(Xs, copy(Xᵒ))
        push!(logweights, logweight)
    end
    Xs, logweights
end

"""
    filter(P::SIRguided, Π, prior, n_ensemble, filter_times; ess_threshold = round(n_ensemble/2; digits=0))

Π: vector of length n_particles containing prior probabilities of initial states
prior: prior on parameter vector in P (a tuple with fieldnames μ, λ and ν)
n_ensemble: number of ensembles in sequential filtering
filter_times: vector of integers that should be strictly increasing, starting with 0 and the final value equal to n_times

Returns n_ensemble paths along with logweight of each path

Suppose P and prior are given then the following gives an example:
n_ensemble = 5
n_particles = nparticles(P)
n_times = ntimes(P)
Π = [SA_F64[0.96, 0.04, 0.0] for _ in 1:n_particles]
Πloc = [Π for _ in 1:n_ensemble]
i0, iT = 3, 5
𝕏, logweights =  guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble, filter_times);
"""
function filter(P::SIRguided, Π, prior, n_ensemble, filter_times; ess_threshold = round(n_ensemble/2; digits=0))
    @assert length(Π)==nparticles(P) "length Π is incorrect"
    @assert filter_times[1]==0 "first element of filter_times should be 0"
    @assert filter_times[end]==ntimes(P) "last element of filter_times should be n_times"

    Πloc = [Π for _ in 1:n_ensemble] # have prior on initial state the same for all ensembles
    i0 = filter_times[1] + 1
    iT = filter_times[2]
    println("Filterig at time $iT")
    𝕏, logweights =  guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble);

    # resampling 
    indices = resamp(logweights, ess_threshold) 
    #println(indices, "\n")
    for i in eachindex(indices)
        𝕏[i] = copy(𝕏[indices[i]])
        logweights[i] = logweights[indices[i]]
    end
    # updating Πloc 
    for i in eachindex(Πloc)
        Πloc[i] = convert_state_to_Π.(𝕏[i][end])  # overwrite element i
    end
    ess_, _= ess(logweights)
    println("ess after resampling equals: ", ess_)

    for i in 3:lastindex(filter_times)
        i0 = filter_times[i-1] + 1
        iT = filter_times[i]
        println("Filterig at time $iT")
        Ys, logweights_ =  guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble);
        println("Interval [$i0, $iT]")
        for j in eachindex(𝕏)
            append!(𝕏[j], Ys[j])
        end
        logweights += logweights_
         # resampling 
        indices = resamp(logweights, ess_threshold) 
        #println(indices, "\n")
        for i in eachindex(indices)
            𝕏[i] = copy(𝕏[indices[i]])
            logweights[i] = logweights[indices[i]]
        end
        # updating Πloc 
        for i in eachindex(Πloc)
            Πloc[i] = convert_state_to_Π.(𝕏[i][end])  # overwrite element i
        end
        ess_, _= ess(logweights)
        println("ess after resampling equals: ", ess_)
    
    end
    𝕏, logweights
end


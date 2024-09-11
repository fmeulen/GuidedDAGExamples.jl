@concrete struct Observation
    t           # time index of observation
    ind         # at each time, indices of observed particles
    h           # corresponding h-vectors
    x           # values at all times (latent = _L_)
end

# Define a custom show method for Observation
function Base.show(io::IO, m::Observation)
    println(io, "Observation Data:")
    println(io, "  t (time index of observation): ", m.t)
    println(io, "  ind (indices of observed particles): ", m.ind)
    println(io, "  h (corresponding h-vectors): ", m.h)
    println(io, "  x (observation, latent = _L_): ", m.x)
end


function observationmessage(x::State)
    if x==_S_
        return(SA_F64[1, 0, 0])
    elseif x==_I_
        return(SA_F64[0, 1, 0])
    elseif x==_R_
        return(SA_F64[0, 0, 1])
    else
        return(SA_F64[1, 1, 1])
    end
end

"""
    create_data(Xtrue, samplesize, n_times, n_particles, O)

    Create a vector of type Observation and length n_times
    Out of n_times * n_particles in total samplesize (time,particle) pairs are randomly chosen
    At those times we observe the process through O
"""
function create_data(Xtrue, samplesize, n_times, n_particles, O)
    # randomly distribute the sample size over times ∈ {1,...,n_times}
    nrobs_at_each_time = rand(Multinomial(samplesize, n_times))  
    ind_obs = Vector{Int64}[]
    for i in 1:n_times
        k = nrobs_at_each_time[i]
        ids = StatsBase.sample(1:n_particles, k; replace=false)   # sample which particles are observed at that time
        push!(ind_obs, ids)
    end

    Xobs = [fill(_L_, n_particles) for _ in 1:n_times   ]
    h = Vector{SVector{3, Float64}}[]
    for t in 1:n_times
        ht = SVector{3, Float64}[]
        for i in ind_obs[t]
            x = Xtrue[t][i]
            #println(x)
            Xobs[t][i] = x
            push!(ht, O * observationmessage(x))
        end
        push!(h, ht)
    end

    [Observation(t, ind_obs[t], h[t], Xobs[t]) for t in 1:n_times]
end

"""
    create_data_regular(Xtrue, obs_times, n_times, n_particles, O)

    Create a vector of type Observation and length n_times
    At all times in obs_times we observe all particles 
    At those times we observe the process through O
"""
function create_data_regular(Xtrue, obs_times, n_times, n_particles, O)
    ind_obs = Vector{Int64}[]
    for i in 1:n_times
        ids = Int64[]
        if i ∈ obs_times
            ids = collect(1:n_particles) # observe all at those times
        end
        push!(ind_obs, ids)
    end

    Xobs = [fill(_L_, n_particles) for _ in 1:n_times   ]
    h = Vector{SVector{3, Float64}}[]
    for t in 1:n_times
        ht = SVector{3, Float64}[]
        for i in ind_obs[t]
            x = Xtrue[t][i]
            Xobs[t][i] = x
            push!(ht, O * observationmessage(x))
        end
        push!(h, ht)
    end

    [Observation(t, ind_obs[t], h[t], Xobs[t]) for t in 1:n_times]
end


"""
    create_data_see_some_particles(Xtrue, obs_particles, n_times, n_particles, O)

    Create a vector of type Observation and length n_times
    At all times in obs_times we observe all particles 
    At those times we observe the process through O
"""
function create_data_see_some_particles(Xtrue, obs_particles, n_times, n_particles, O)
    ind_obs = fill(collect(obs_particles), n_times)
    
    Xobs = [fill(_L_, n_particles) for _ in 1:n_times   ]
    h = Vector{SVector{3, Float64}}[]
    for t in 1:n_times
        ht = SVector{3, Float64}[]
        for i in ind_obs[t]
            x = Xtrue[t][i]
            Xobs[t][i] = x
            push!(ht, O * observationmessage(x))
        end
        push!(h, ht)
    end

    [Observation(t, ind_obs[t], h[t], Xobs[t]) for t in 1:n_times]
end

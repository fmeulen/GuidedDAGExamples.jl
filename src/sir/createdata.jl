@concrete struct Observation
    ind         # at each time, indices of observed particles
    h           # corresponding h-vectors
    x           # values at all times (latent = _L_)
end

function create_data(samplesize, n_times, n_particles, O)
    nrobs_at_each_time = rand(Multinomial(samplesize, n_times))  
    ind_obs = Vector{Int64}[]
    for i in 1:n_times
        k = nrobs_at_each_time[i]
        ids = StatsBase.sample(1:n_particles, k; replace=false)   
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

    [Observation(ind_obs[t], h[t], Xobs[t]) for t in 1:n_times]
end


# function show(::Observation, 𝒪)
#     for x in 𝒪 
#         println(x.ind)
#     end
# end
function construct_Ploc(P, i0, iT)
    𝒪loc = P.𝒪[i0:iT]
#    ℐloc = P.ℐ[i0:iT]   
    ℐloc =initialise_infected_neighbours(𝒪loc) # goes wrong if there are no observations in time slice
    SIRguided(P.ξ, P.λ, P.μ, P.ν, P.τ, 𝒩, ℐloc, 𝒪loc, O) 
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

function guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble)
    Ploc = construct_Ploc(P,i0,iT)
    Bloc = backward(Ploc)

    n_times = iT - i0 + 1
    n_particles = length(P.𝒩)
    Z = innovations(n_times, n_particles)
    Xᵒ, logweight  = forward(Ploc, Πloc[1], Bloc, Z, prior)

    Xs = [Xᵒ]
    logweights = [logweight]
    
    # repeat this for all particles
    for i in 2:n_ensemble 
        Z = innovations(n_times, n_particles)
        Xᵒ, logweight  = forward(Ploc, Πloc[i], Bloc, Z, prior)
        push!(Xs, Xᵒ)
        push!(logweights, logweight)
    end
    println("logweights: ", logweights)

    # resampling
    ess_threshold = round(n_ensemble/2; digits=0) 
    indices = resamp(logweights, ess_threshold) 
    println(indices, "\n")
    for i in eachindex(indices)
        Xs[i] = Xs[indices[i]]
        logweights[i] = logweights[indices[i]]
    end
    
    for i in eachindex(Πloc)
        Πloc[i] = convert_state_to_Π.(Xs[i][end])  # overwrite element i
    end
    Xs, logweights
end


lo = @layout [a;b;c]

n_ensemble = 20
filter_times = vcat(0,20:10:100)
Πloc = [Π for _ in 1:n_ensemble]
i = 2
i0 = filter_times[i-1] + 1
iT = filter_times[i]
Xs, logweights =  guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble);

for i in 3:length(filter_times)
    i0 = filter_times[i-1] + 1
    iT = filter_times[i]
    println("i0, iT= $i0, $iT")
    Ys, logweights =  guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble);
    for i in eachindex(Xs)
        append!(Xs[i], Ys[i])
    end
end

plotpath(Xs[1]; xlims_=(1,n_times))

plot(pforward, plotpath(Xs[1]), pobs, layout=lo)

plotpath(particles[1].X)
xlims=(1,n_times)

lo = @layout [a;b;c]
plot(pforward,plotpath(particles[1].X;name="guided"),pobs, layout=lo)

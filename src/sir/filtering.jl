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
        push!(Xs, copy(Xᵒ))
        push!(logweights, logweight)
    end
 #   println("logweights: ", logweights)

    # resampling
    ess_threshold = round(n_ensemble/2; digits=0) 
    indices = resamp(logweights, ess_threshold) 
    println(indices, "\n")
    for i in eachindex(indices)
        Xs[i] = copy(Xs[indices[i]])
        logweights[i] = logweights[indices[i]]
    end
    
    for i in eachindex(Πloc)
        Πloc[i] = convert_state_to_Π.(Xs[i][end])  # overwrite element i
    end
    Xs, logweights
end


lo = @layout [a;b;c]



# Suppose P and prior are given then the following gives an example:
n_ensemble = 5
n_particles = nparticles(P)
n_times = ntimes(P)
Π = [SA_F64[0.96, 0.04, 0.0] for _ in 1:n_particles]
Πloc = [Π for _ in 1:n_ensemble]
i0, iT = 3, 5
Xs, logweights =  guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble);

for j in 1:length(Xs)
    println("Xs[1] === Xs[$j]: ", Xs[1] === Xs[j])
end

function filterall(P, Π, prior, n_ensemble)
    Πloc = [Π for _ in 1:n_ensemble]

    i = 2
    i0 = filter_times[i-1] + 1
    iT = filter_times[i]
    𝕏, logweights =  guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble);

    for i in 3:lastindex(filter_times)
        println("i=$i, length XXs[1]=", length(XXs[1]))
        i0 = filter_times[i-1] + 1
        iT = filter_times[i]
        println("i0, iT= $i0, $iT")
        Ys, logweights =  guidedfiltering!(Πloc, P, i0, iT, prior, n_ensemble);
        println("Interval [$i0, $iT], Length of Ys[1]: ", length(Ys[1]))
        for j in eachindex(𝕏)
            append!(𝕏[j], Ys[j])
            println("length of 𝕏[j] ", length(𝕏[j]))
            #push!(XXs, Ys)
        end
    end
    𝕏
end

n_ensemble = 100
filter_times = [0, 30, 50, 100]

𝕏 = filterall(P, Π, prior, n_ensemble);
[length(𝕏[i]) for i in eachindex(𝕏)]

anim_smc = @animate for  x ∈ 𝕏
    #plotpath(𝕏[1]; xlims_=(1,n_times))
    plot(pforward, plotpath(x), pobs, layout=lo)
end

mp4(anim_smc,presfigdir*"/smc_example.mp4", fps=1.5)

xlims=(1,n_times)

lo = @layout [a;b;c]
plot(pforward,plotpath(particles[1].X;name="guided"),pobs, layout=lo)

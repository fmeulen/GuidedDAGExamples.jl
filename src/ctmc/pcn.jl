
z = randn(1000)
rexp(z) = -logcdf.(Normal(), z)
x = rexp(z)

R = RNᵒ

function sample_pcn(t0,x0,z, R::Network; Tend=nothing)
    t = t0
    x = x0
    T = zeros(R.nreact)
    counter = fill(2,R.nreact)

    P = rexp([z[i][1] for i ∈ 1:R.nreact])
    eventtimes = [t]
    eventvals = [x]
    if isnothing(Tend)  Tend=R.T  end
    while t < Tend
        T, P, x, t, z, counter = simstep_pcn!(T, P, x, t, z, counter, R)
        push!(eventtimes, t)
        push!(eventvals, x)
    end
    if last(eventtimes) >= Tend
        deleteat!(eventtimes,length(eventtimes))
        deleteat!(eventvals,length(eventvals))
    end
    eventtimes, eventvals, z
end

function simstep_pcn!(T, P, x, t, z, counter, R::GuidedReactionNetwork; ρ=0.7)
    aᵒ = λᵒ(t,x,R)
    Δt = abs.((P - T) ./ aᵒ)
    (Δ, μ) = findmin(Δt)
    x += R.ξ[μ]   # update state according to reaction taking place
    t += Δ  # update clock time
    T += Δ * aᵒ  # update internal times
    zᵒ = ρ * z[μ][counter[μ]] + sqrt(1-ρ^2) * randn()  # pCN step
    P[μ] += rexp(zᵒ)
    counter[μ] += 1
    T, P, x, t, z, counter
end

RN = ReactionNetwork(1.5, 3.0, [[1,1], [-1,0], [0,1], [0,-1]],4)
Tend = 53.0
x0 = [16,5]

times_forw, events_forw = sample(t0,x0, RN;Tend=Tend)
xT = last(events_forw)


μ = 0.0 * μT(Tend, x0, xT)
Σ = ΣT(RN,xT)
RNᵒ = guidedreactionnetwork(RN,Tend,xT,μ,Σ)
times_guid, events_guid = sample(t0,x0, RNᵒ)


z = [randn(1000000) for _ ∈ 1:RNᵒ.nreact]
times_guid, events_guid, z = sample_pcn(t0,x0,z, RNᵒ)
ll = loglik(times_guid, events_guid, RNᵒ)

function mh_pcn(t0, x0, RNᵒ, iter)
    z = [randn(10000) for _ ∈ 1:RNᵒ.nreact]
    times_guid, events_guid, z = sample_pcn(t0,x0,z, RNᵒ)
    ll = loglik(times_guid, events_guid, RNᵒ)

    Z = [z]
    LL = [ll]
    T = [times_guid]
    E = [events_guid]

    acc = 0

    for j in 1:iter
        times_guidᵒ, events_guidᵒ, zᵒ = sample_pcn(t0,x0,z, RNᵒ)
        llᵒ = loglik(times_guidᵒ, events_guidᵒ, RNᵒ)
        println(llᵒ-ll)
        if log(rand()) < llᵒ - ll
            ll = llᵒ
            z = zᵒ
            push!(T, times_guidᵒ)
            push!(E, events_guidᵒ)
            acc += 1

        else
            push!(T, times_guid)
            push!(E, events_guid)
        end
        push!(LL,ll)
    end
    println("Average acceptance: $acc")
    T, E, LL
end

T, E, LL = mh_pcn(t0, x0, RNᵒ, 1000)



p = plot(times_forw, first.(events_forw),label="forward el1", legend = :outertopleft)
plot!(p, times_forw, last.(events_forw),label="forward el2")

#for i ∈ Set([2,4,10])
for i ∈ [2, 1000]
    plot!(p, T[i], first.(E[i]),label="guided el1")
    plot!(p, T[i], last.(E[i]),label="guided el2")
end
p
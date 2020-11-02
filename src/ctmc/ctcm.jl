using Distributions
using Plots
using Random
using LinearAlgebra
using Roots
using RCall
using SpecialFunctions
#@rlibrary pracma

include("rootsolving.jl")

abstract type Network end

struct ReactionNetwork <: Network
    λ::Float64
    ρ::Float64
    ξ::Array{Array{Int64,1},1} # transitions
    nreact::Int64
end

struct GuidedReactionNetwork <: Network
    λ::Float64
    ρ::Float64
    ξ::Array{Array{Int64,1},1}
    nreact::Int64
    T::Float64 # endtime
    xT::Array{Int64,1} # endpoint
    μ::Array{Float64,1}  # mean of X̃
    Σ::Array{Float64,2} # cov matrix of X̃
end

guidedreactionnetwork(RN::ReactionNetwork,T,xT,μ,Σ) = GuidedReactionNetwork(RN.λ, RN.ρ, RN.ξ, RN.nreact, T, xT, μ, Σ)

"""
    simstep!(T, P, x, t, R::ReactionNetwork)

one step of next-reaction method
"""
function simstep!(T, P, x, t, R::ReactionNetwork)
    a = λ(x,R)
    Δt = abs.((P - T) ./ a)
    (Δ, μ) = findmin(Δt)
    t += Δ  # update clock time
    T += Δ * a  # update internal times
    x += R.ξ[μ]   # update state according to reaction taking place
    P[μ] += randexp()
    T, P, x, t
end


"""
    sample(t0, Tend, x0, R::Network)

Sample network starting at time t0 in state x0 up to time Tend
Returns eventtimes and eventvalues
"""
function sample(t0, Tend, x0, R::Network)
    t = t0
    x = x0
    T = Base.zeros(R.nreact)
    P = randexp(R.nreact)
    eventtimes = [t]
    eventvals = [x]
    while t < Tend
        T, P, x, t = simstep!(T, P, x, t, R)
        push!(eventtimes, t)
        push!(eventvals, x)
    end
    if last(eventtimes) >= Tend
        deleteat!(eventtimes,length(eventtimes))
        deleteat!(eventvals,length(eventvals))
#        push!(eventtimes, Tend)
#        push!(eventvals, last(eventvals))  # can also choose to put xT here
    end
    eventtimes, eventvals
end

#####################################################################
# Example forward simulation

RN = ReactionNetwork(0.5, 6.0, [[1,0], [-1,0], [0,1], [0,-1]],4)
λ(x, RN::Network) =[RN.ρ/(1.0+x[2]), RN.λ*(x[1] >= 1), RN.ρ/(1.0+x[1]), RN.λ*(x[2] >= 1)]

t0 = 0.0
x0 = [2, 0]
Tend =30
eventtimes, eventvals = sample(t0, Tend, x0, RN)

p = plot(eventtimes, first.(eventvals))
plot!(p, eventtimes, last.(eventvals))

#####################################################################
# Guiding
ψ(s, x, R::GuidedReactionNetwork) = [LinearAlgebra.dot(-0.5R.ξ[ℓ] + R.xT - (x + R.μ *(R.T - s)), R.Σ\R.ξ[ℓ]) for ℓ ∈ eachindex(R.ξ)]  # assumigng Σ is identity matrix
λᵒ(s, x, R::GuidedReactionNetwork) = λ(x, R) .* exp.(ψ(s,x,R)/(R.T - s))
logh̃(s,x,R::GuidedReactionNetwork) = logpdf(MvNormal(x + R.μ*(R.T-s), R.Σ*(R.T - s)), R.xT)

"""
     simstep!(T, P, x, t, R::GuidedReactionNetwork; approximate=true)

One step of next-reaction method for guided reaction network.

Approximate method treats the time-dependent density of time-independent, by plugging in the left-value of t
The exact method is still buggy.
"""
function simstep!(T, P, x, t, R::GuidedReactionNetwork; approximate=true) #FIXME
    aᵒ = λᵒ(t,x,R)
    if approximate
        Δt = abs.((P - T) ./ aᵒ)
    else
        Ψ = ψ(t,x,R)
        L = λ(x,R)
        Δt = [rootsolve(L[ℓ], R.T - t, Ψ[ℓ], P[ℓ]-T[ℓ]) for ℓ ∈ 1:R.nreact]
    end
    (Δ, μ) = findmin(Δt)
    x += R.ξ[μ]   # update state according to reaction taking place
    t += Δ  # update clock time
    if approximate
        Ψ = ψ(0.0, x, R)
        T += λ(x,R) .* [g(R.T - t, R.T - t + Δ, -Ψ[ℓ]) for ℓ ∈ 1:R.nreact]  # update internal times
    else
        T += Δ * aᵒ
    end
    P[μ] += randexp()
    T, P, x, t
end

""""
    loglik(t,v, Rᵒ)

(t,v) = (eventtimes, eventvals)
"""
function loglik(t,v, Rᵒ; sk=1) #FIXME
    ll = logh̃(t[1],v[1],Rᵒ)
    for k ∈ 1:length(t)-1-sk
        ll += logh̃(t[k+1],v[k],Rᵒ) - logh̃(t[k],v[k],Rᵒ)
        #ll += (t[k+1]-t[k]) * sum(λᵒ(t[k],v[k], Rᵒ) - λ(v[k], Rᵒ))
        L = λ(v[k], Rᵒ)
        Ψ = ψ(0.0, v[k], Rᵒ)
        for ℓ ∈ 1:Rᵒ.nreact
            if iszero(Ψ[ℓ])
                intcontrib = 0.0
            else
                intcontrib = g(Rᵒ.T-t[k+1], Rᵒ.T-t[k], -Ψ[ℓ]) - (t[k+1]-t[k])
            end
            if isnan(intcontrib)
                println("k=$k, ℓ=$ℓ, Ψ=$Ψ")
            end
            Ψℓ = Ψ[ℓ]
            # problems if Ψℓ>0 on very last interval
        #    println("k=$k, ℓ=$ℓ, intcontrib=$intcontrib, Ψℓ = $Ψℓ")
            ll += L[ℓ] * intcontrib
        end
    end
    ll
end

#####################################################################

RN = ReactionNetwork(1.5, 3.0, [[1,1], [-1,0], [0,1], [0,-1]],4)
Tend = 33.0
x0 = [16,5]

function ΣT(R::ReactionNetwork, xT)
    LT = λ(xT, R)
    Σ = zeros(length(xT), length(xT))
    for ℓ ∈ eachindex(R.ξ)
        Σ += R.ξ[ℓ] * R.ξ[ℓ]' * LT[ℓ]
    end
    Σ
end

μT(Tend, x0, xT) = (xT-x0)/Tend

times_forw, events_forw = sample(t0, Tend, x0, RN)
xT = last(events_forw)

μ = 0.0 * μT(Tend, x0, xT)
Σ = .1 * ΣT(RN,xT)  # this choice of open to discussion (a smaller multiplicative factor
# (here 0.1) causes a stronger pull, and, once the endpoint is reached, to stay there)

h = 0.001 # added to shift conditioning time slightly beyond T
RNᵒ = guidedreactionnetwork(RN, Tend  +h, xT ,μ, Σ)
times_guid, events_guid = sample(t0, Tend, x0, RNᵒ)
ll = loglik(times_guid, events_guid, RNᵒ)

print(last(events_guid), xT)

p = plot(times_forw, first.(events_forw),label="forward el1", legend = :outertopleft)
plot!(p, times_forw, last.(events_forw),label="forward el2")
plot!(p, times_guid, first.(events_guid),label="guided el1")
plot!(p, times_guid, last.(events_guid),label="guided el2")

ll = loglik(times_forw, events_forw, RNᵒ)
ll = loglik(times_guid, events_guid, RNᵒ)

# get guiding intensities of guided proposal at all event times
#[λ(times_guid[i], events_guid[i], RNᵒ) for i ∈ eachindex(times_guid)]

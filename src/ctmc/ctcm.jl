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
    xT::Array{Float64,1} # endpoint # FIXME depending on noise Vector{Int64}
    μ::Array{Float64,1}
    Σ::Array{Float64,2} # noise in guiding term
end

guidedreactionnetwork(RN::ReactionNetwork,T,xT,μ,Σ) = GuidedReactionNetwork(RN.λ, RN.ρ, RN.ξ, RN.nreact, T, xT, μ, Σ)

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
        push!(eventtimes, Tend)
        push!(eventvals, last(eventvals))  # can also choose to put xT here
    end
    eventtimes, eventvals
end

#####################################################################
# Forward simulation

# specify reactions
#RN = ReactionNetwork(0.5, 6.0, [[1,0], [-1,0], [0,1], [0,-1]],4)
RN = ReactionNetwork(0.5, 6.0, [[2,1], [-1,0], [0,1], [0,-1]],4)
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

function simstep!(T, P, x, t, R::GuidedReactionNetwork; approximate=true)#false)
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
    Ψ = ψ(0.0, x, R)
    #T += λ(x,R) .* [g(R.T - t, R.T - t + Δ, -Ψ[ℓ]) for ℓ ∈ 1:R.nreact]  # update internal times
    T += Δ * aᵒ
    P[μ] += randexp()
    T, P, x, t
end

""""
    loglik(t,v, Rᵒ)

(t,v) = (eventtimes, eventvals)
"""
function loglik(t,v, Rᵒ) #FIXME
    ll = logh̃(t[1],v[1],Rᵒ)
    for k ∈ 1:length(t)-1
        #if t[k+1] >= Rᵒ.T @error "eventtime larger than T"  end
        ll += logh̃(t[k+1],v[k],Rᵒ) - logh̃(t[k],v[k],Rᵒ)
        ll += (t[k+1]-t[k]) * sum(λᵒ(t[k],v[k], Rᵒ) - λ(v[k], Rᵒ))
    end
    ll
end
#####################################################################

RN = ReactionNetwork(1.5, 3.0, [[1,1], [-1,0], [0,1], [0,-1]],4)
Tend = 33.0
x0 = [16,5]
xT = [20,3]

Σ = 0.0 *  Matrix(1.0I,2,2)

function ΣT(R::ReactionNetwork, xT)
    LT = λ(xT, R)
    Σ = zeros(length(xT), length(xT))
    for ℓ ∈ eachindex(R.ξ)
        Σ += R.ξ[ℓ] * R.ξ[ℓ]' * LT[ℓ]
    end
    Σ
end

μT(Tend, x0, xT) = (xT-x0)/Tend

Tend = 40.0
times_forw, events_forw = sample(t0,x0, RN;Tend=Tend)
xT = last(events_forw) #+ 0.001*fill(1.0,2)#*randn(length(x0))

μ = 0.0 * μT(Tend, x0, xT)
Σ = 0.1 * ΣT(RN,xT)
RNᵒ = guidedreactionnetwork(RN,Tend,xT,μ,Σ)
times_guid, events_guid = sample(t0, Tend, x0, RNᵒ)
#ll = loglik(times_guid, events_guid, RNᵒ)

print(last(events_guid), xT)

p = plot(times_forw, first.(events_forw),label="forward el1", legend = :outertopleft)
plot!(p, times_forw, last.(events_forw),label="forward el2")
plot!(p, times_guid, first.(events_guid),label="guided el1")
plot!(p, times_guid, last.(events_guid),label="guided el2")



# get guiding intensities of guided proposal at all event times
#[λ(times_guid[i], events_guid[i], RNᵒ) for i ∈ eachindex(times_guid)]


if false

# new try
#ζ(x, R::GuidedReactionNetwork) = [LinearAlgebra.dot(R.xT-x,R.ξ[ℓ]) for ℓ ∈ eachindex(R.ξ)]
ζ(s, x, R::GuidedReactionNetwork) = [LinearAlgebra.dot(R.xT - x - (R.T-s)*R.μ, R.ξ[ℓ]) /λ(R.xT, R)[ℓ] for ℓ ∈ eachindex(R.ξ)]
λ(s, x, R::GuidedReactionNetwork) = λ(x, R) .* ζ(s, x, R) /(R.T - s)

function simstep!(T, P, x, t, R::GuidedReactionNetwork)
    z = ζ(t, x, R)
    imin = findall(x-> x < 0, z)
    ipos = findall(x-> x >= 0, z)
    a = λ(x,R)
    #aᵒ = λ(t,x,R)
    Δt = abs.((P - T) ./ a)
    for ℓ ∈ ipos
        Δt[ℓ] = (R.T-t) * (1-exp(-Δt[ℓ]/z[ℓ]))
    end
    (Δ, μ) = findmin(Δt)
    for ℓ ∈ imin
        T[ℓ] += Δ * a[ℓ]  # update internal times ## ADJUST
    end
    for ℓ ∈ ipos
        T[ℓ] += a[ℓ] * z[ℓ] *log((R.T-t)/(R.T-t-Δ))
    end
    t += Δ  # update clock time
    x += R.ξ[μ]   # update state according to reaction taking place
    P[μ] += randexp()
    T, P, x, t
end

end
function τleap(t0,x0, R::GuidedReactionNetwork; h=0.1)
    t = t0
    x = x0
    eventtimes = [t]
    eventvals = [x]
    while t < R.T
        L = λ(t,x,R)
    #    println(L)
        for ℓ ∈ 1:R.nreact
            x += R.ξ[ℓ] * Random.rand(Poisson(h*L[ℓ]))
        end
        t += h
        push!(eventtimes, t)
        push!(eventvals, x)
    end
    eventtimes, eventvals
end

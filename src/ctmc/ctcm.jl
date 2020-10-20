using Distributions
using Plots
using Random
using LinearAlgebra

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
    Σ::Array{Float64,2} # noise in guiding term
end

guidedreactionnetwork(RN::ReactionNetwork,T,xT,Σ) = GuidedReactionNetwork(RN.λ, RN.ρ, RN.ξ, RN.nreact, T, xT, Σ)

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

function simsteps(t0,x0, R::Network; Tend=nothing)
    t = t0
    x = x0
    T = zeros(R.nreact)
    P = randexp(R.nreact)
    eventtimes = [t]
    eventvals = [x]
    if isnothing(Tend)  Tend=R.T  end
    while t < Tend
        T, P, x, t = simstep!(T, P, x, t, R)
        push!(eventtimes, t)
        push!(eventvals, x)
    end
    eventtimes, eventvals
end

#####################################################################
# Forward simulation

# specify reactions
RN = ReactionNetwork(0.5, 6.0, [[1,0], [-1,0], [0,1], [0,-1]],4)
λ(x, RN::Network) =[RN.ρ/(1.0+x[2]), RN.λ*(x[1] >= 1), RN.ρ/(1.0+x[1]), RN.λ*(x[2] >= 1)]

t0 = 0.0
x0 = [2, 0]
eventtimes, eventvals = simsteps(t0,x0, RN; Tend=30)

p = plot(eventtimes, first.(eventvals))
plot!(p, eventtimes, last.(eventvals))


#####################################################################
# Guiding

h̃(s,x,R::GuidedReactionNetwork) = logpdf(MvNormal(x, R.Σ*(R.T-s)), R.xT)
λ(t, x, R::GuidedReactionNetwork) = λ(x, R) .* [h̃(t,x+R.ξ[ℓ],R) for ℓ ∈ eachindex(R.ξ)] /h̃(t,x,R)
ψ(x, R::GuidedReactionNetwork) = [dot(0.5R.ξ[ℓ] - R.xT + x, R.Σ\R.ξ[ℓ]) for ℓ ∈ eachindex(R.ξ)]  # assumigng Σ is identity matrix

function simstep!(T, P, x, t, RNᵒ::GuidedReactionNetwork)
    aᵒ = λ(t,x,RNᵒ)
    a = λ(x,RNᵒ) .* exp.(-ψ(x, RNᵒ)/(Tend - t))
    Δt = abs.((P - T) ./ a)
    (Δ, μ) = findmin(Δt)
    t += Δ  # update clock time
    T += Δ * aᵒ  # update internal times
    x += RNᵒ.ξ[μ]   # update state according to reaction taking place
    P[μ] += randexp()
    T, P, x, t
end

#####################################################################

Tend = 10.0
xT = [15,8]
Σ = 2.0*Matrix(1.0I,2,2)
RNᵒ = guidedreactionnetwork(RN,Tend,xT,Σ)

eventtimes, eventvals = simsteps(t0,x0, RNᵒ)
print(last(eventvals), xT)

p = plot(eventtimes, first.(eventvals))
plot!(p, eventtimes, last.(eventvals))

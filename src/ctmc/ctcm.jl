using Distributions
using Plots
using Random
using LinearAlgebra
using Roots
using RCall
using SpecialFunctions
#@rlibrary pracma

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

function sample(t0,x0, R::Network; Tend=nothing)
    t = t0
    x = x0
    T = Base.zeros(R.nreact)
    P = randexp(R.nreact)
    eventtimes = [t]
    eventvals = [x]
    if isnothing(Tend)  Tend=R.T  end
    while t < Tend
        T, P, x, t = simstep!(T, P, x, t, R)
        push!(eventtimes, t)
        push!(eventvals, x)
    end
    if last(eventtimes) >= Tend
        deleteat!(eventtimes,length(eventtimes))
        deleteat!(eventvals,length(eventvals))
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
eventtimes, eventvals = sample(t0,x0, RN; Tend=30)

p = plot(eventtimes, first.(eventvals))
plot!(p, eventtimes, last.(eventvals))


#####################################################################
# Guiding
ψ(s, x, R::GuidedReactionNetwork) = [LinearAlgebra.dot(-0.5R.ξ[ℓ] + R.xT - (x + R.μ *(R.T - s)), R.Σ\R.ξ[ℓ]) for ℓ ∈ eachindex(R.ξ)]  # assumigng Σ is identity matrix
λᵒ(s, x, R::GuidedReactionNetwork) = λ(x, R) .* exp.(ψ(s, x,R)/(R.T-s))
logh̃(s,x,R::GuidedReactionNetwork) = logpdf(MvNormal(x + R.μ*(R.T-s), R.Σ*(R.T - s)), R.xT)

#--------- begin Moritz' stuff
Ei(x) = -real(expint(Complex(-x)))
g(u, a) = exp(-a/u)*u + a*Ei(-a/u)

"""
    g(l, u, a)

Integrate  exp(-a/x) from l to u
"""
g(l, u, a) = g(u, a) - g(l, a)
g(1, 2, pi) # should be 0.124
# https://www.wolframalpha.com/input/?i=integrate+exp%28-pi%2Fu%29+from+1+to+2
#---------end  Moritz' stuff


# ei(x) = expint(x)[1]
# ei1(x)= expint_E1(x)[1]
# g(ℓ, u, α) = α>0 ? u*exp(-α/u) - ℓ*exp(-α/ℓ) + α*(ei(α/u) - ei(α/ℓ)) : u*exp(-α/u) - ℓ*exp(-α/ℓ) - α*(ei1(-α/ℓ) - ei1(-α/u))
G(λℓ, ψℓ, Tmint, PℓminTℓ) = (Δ) -> λℓ * g(Tmint-Δ, Tmint, -ψℓ) - PℓminTℓ
∇G(λℓ, ψℓ, t,R::GuidedReactionNetwork) = (Δ) -> -λℓ *   exp(ψℓ/(R.T-t-Δ))

function simstep!(T, P, x, t, R::GuidedReactionNetwork; approximate=true)#false)
    aᵒ = λᵒ(t,x,R)
    if approximate
        Δt = abs.((P - T) ./ aᵒ)
    else
        Δt = zeros(R.nreact)
        L = λ(x,R)
        Ψ = ψ(t,x,R)
        for ℓ ∈ eachindex(Δt)
            ff = (G(L[ℓ], Ψ[ℓ], R.T-t, P[ℓ]-T[ℓ]), ∇G(L[ℓ],Ψ[ℓ] ,t, R))
            println((L[ℓ], Ψ[ℓ], R.T-t, P[ℓ]-T[ℓ]))
            Δt[ℓ] = find_zero(ff, 0.0001,  Roots.Newton())      # (0.0001,R.T-t); maxevals =50)#,R.T-t) # returns Float64
        end
    end
    (Δ, μ) = findmin(Δt)
    x += R.ξ[μ]   # update state according to reaction taking place
    t += Δ  # update clock time
    T += Δ * aᵒ  # update internal times
    P[μ] += randexp()
    T, P, x, t
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

""""
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
Tend = 13.0
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


times_forw, events_forw = sample(t0,x0, RN;Tend=Tend)
xT = last(events_forw)

μ = 0.0 * μT(Tend, x0, xT)
Σ = ΣT(RN,xT)#*Tend/sum(xT))#Σ =  1*Matrix(1.0I,2,2)
RNᵒ = guidedreactionnetwork(RN,Tend,xT,μ,Σ)
times_guid, events_guid = sample(t0,x0, RNᵒ)
ll = loglik(times_guid, events_guid, RNᵒ)

print(last(events_guid), xT)

p = plot(times_forw, first.(events_forw),label="forward el1", legend = :outertopleft)
plot!(p, times_forw, last.(events_forw),label="forward el2")
plot!(p, times_guid, first.(events_guid),label="guided el1")
plot!(p, times_guid, last.(events_guid),label="guided el2")


# Independence MH
res = [[times_guid, events_guid, ll]]
for i=1:1000
    global ll, times_guid, events_guid, res
    times_guidᵒ, events_guidᵒ = sample(t0,x0, RNᵒ)
    llᵒ = loglik(times_guidᵒ, events_guidᵒ, RNᵒ)
#    println(llᵒ-ll)
    if log(rand()) < llᵒ - ll
                push!(res,[times_guidᵒ, events_guidᵒ, llᵒ])
        ll = llᵒ
    else
        push!(res,[times_guid, events_guid, ll])
    end
end

times_guid, events_guid = last(res)[1], last(res)[2]
print(last.(res))

p = plot(times_forw, first.(events_forw),label="forward el1", legend = :outertopleft)
plot!(p, times_forw, last.(events_forw),label="forward el2")
plot!(p, times_guid, first.(events_guid),label="guided el1")
plot!(p, times_guid, last.(events_guid),label="guided el2")


println(times_guid[end-5:end])
println(eventvals[end-5:end])

# get guiding intensities of guided proposal at all event times
[λ(times_guid[i], events_guid[i], RNᵒ) for i ∈ eachindex(times_guid)]


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

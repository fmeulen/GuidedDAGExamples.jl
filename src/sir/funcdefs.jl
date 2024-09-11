using StaticArrays
abstract type MarkovProcess end

# abstract type BackwardFilter end
# struct Diagonal <:


@concrete struct SIRforward <: MarkovProcess
        ξ::Float64  # prob to stay in susceptible state if there are no infected neighbours
        λ::Float64
        μ::Float64
        ν::Float64
        τ::Float64
        𝒩::Array{Array{Int64,1},1}
end

@concrete struct SIRguided{T, S} <: MarkovProcess
    ξ::Float64
    λ::Float64
    μ::Float64
    ν::Float64
    τ::Float64
    𝒩::Array{Array{Int64,1},1}  # neighbourhood structure
    ℐ::Array{Array{T,1},1} # vector with each element a vector at that time of nr of infected neighbours
    𝒪::Vector{S}
    O::SMatrix{3, 3, Float64, 9} 
end

SIRguided(P::SIRforward, ℐ, 𝒪, O) = SIRguided(P.ξ, P.λ , P.μ, P.ν, P.τ, P.𝒩, ℐ, 𝒪, O)

param(P::SIRguided) = (λ=P.λ, μ=P.μ, ν=P.ν)

@enum State::UInt8 _S_=1 _I_=2 _R_=3 _L_=0
const 𝒳 = @SVector [_S_,_I_,_R_]

# state space of one person
const BF = [@SVector([3,1,2]), @SVector([1,2,3]), @SVector([2,3,1])]


######## function defns

params(P::MarkovProcess) = (λ=P.λ, μ=P.μ, ν=P.ν)

"""
    ind(x)

x ∈ 𝒳 is mapped into integer
"""
ind(x) = Int(x)



"""
    set_neighbours(n)

Returns an array of length n, where the i-th element contains the indices for the neighbours of individual i
"""
function set_neighbours(n, k)
    if k==2
        if n<5
        @error "provide larger value for n"
        end
        𝒩 = Vector{Vector{Int64}}(undef,n)
        𝒩[1] = [2,3]; 𝒩[2] = [1,3,4]; 𝒩[n-1]  = [n, n-2, n-3]; 𝒩[n] = [n-1,n-2]
        for i in 3:n-2
            𝒩[i] = [i-2,i-1,i+1,i+2]
        end
    end
    if k==1
        𝒩 = Vector{Vector{Int64}}(undef,n)
        𝒩[1] = [2]; 𝒩[n] =[n-1]
        for i in 2:n-1
            𝒩[i] = [i-1, i+1]
        end
    end
    𝒩
end

#### forward sampling #####
rand𝒳(z, p)= z < p[1] ? _S_ : ( z > 1-p[3] ? _R_ : _I_ )

function rand𝒳!(x::State, z, p)
    if z < p[1]
        x= _S_
    elseif z > 1-p[3]
        x=  _R_
    else
        x= _I_
    end
end

# """
#     sample𝒳(x,z,p)

# Sample from 𝒳 according to probability vector p

# p =  [1/4, 1/2, 1/4]
# z = rand()
# sample𝒳(_S_, z, p)
# """
# function sample𝒳(x, u::Float64,p) # provide current state
#     b, f = BF[Int(x)][1], BF[Int(x)][3]
#     if u < p[b] return(𝒳[b])
#     elseif u > 1-p[f] return(𝒳[f])
#     else return x
#     end
# end

# # old 
# pS(x) = @SVector [exp(-x), 1.0-exp(-x), 0.0]  # λ*τ
# pI(x) = @SVector [0.0, exp(-x), 1.0-exp(-x)]      # μ*τ
# pR(x) = @SVector [1.0-exp(-x), 0.0, exp(-x)]     # ν*τ



# function κ(P, i, x) # used to be called pi
#     if x[i] == _S_
#         pS(P.λ * nr_infected_neighb(x,P.𝒩,i) * P.τ)
#     elseif x[i] == _I_
#         pI(P.μ * P.τ)
#     elseif x[i] == _R_
#         pR(P.ν * P.τ)
#     end
# end

# old 
function pS(P, ninfected) 
    y = P.ξ * exp(-P.λ * ninfected * P.τ)
    @SVector [y,  1.0-y, 0.0]  
end

function pI(P) 
    y = exp(-P.μ * P.τ)
    @SVector [0.0, y, 1.0-y]      # μ*τ
end 


function pR(P) 
    y = exp(-P.ν * P.τ)
    @SVector [1.0-y, 0.0, y]     # ν*τ
end

# returns probability vector for i-th particle if current state (of all particles) is x
function κ(P, i, x) 
    if x[i] == _S_
        ninfected = nr_infected_neighb(x, P.𝒩, i)
        return(pS(P, ninfected))
    elseif x[i] == _I_
        return(pI(P))
    elseif x[i] == _R_
        return(pR(P))
    end
end

# backward kernel for one individual (a matrix)
κ̃(P ,ninfected::Number) = hcat(pS(P, ninfected), pI(P), pR(P))'



"""
    sample_particle(P::SIRforward,i,x,z)

One forward simulation step for invdividual `i`, if the present configuration of all individuals is `x`
i: index of individual to be forward simulated
x: current state of all individuals
z: innovation for this step
"""
function sample_particle(P::SIRforward, i, x, z)
    p = κ(P, i, x)
    rand𝒳(z, p)  #    sample𝒳(x[i], z, p)
end

sample(P::SIRforward, x, z) = [sample_particle(P, i, x, z[i]) for i in eachindex(x)]

"""
    sample_trajectory(P::SIRforward, n::Int, x0)

    sample SIR-process over n time instances, where the first time instance is x0
"""
function sample_trajectory(P::SIRforward, n_times::Int64, x1) 
    X = [x1]
    n_particles = length(x1)
    for j in 2:n_times
        z = rand(n_particles)
        push!(X, sample(P, X[j-1], z))
    end
    X
end







obs2matrix(X) =  [ind(X[j][i]) for j in eachindex(X), i in eachindex(X[1])]


# Define a custom show method for SIRforward
function Base.show(io::IO, m::SIRforward)
    println(io, "SIRforward Process:")
    println(io, "  ξ (1-ξ is probability of new exogeneous infection): ", m.ξ)
    println(io, "  λ (infection rate): ", m.λ)
    println(io, "  μ (recovery rate): ", m.μ)
    println(io, "  ν (rate to get susceptible again): ", m.ν)
    println(io, "  τ (discretisation time step): ", m.τ)
    println(io, "  𝒩 (neighbors list): ")
    for (i, neighbors) in enumerate(m.𝒩)
        println(io, "    Neighbor $i: ", neighbors)
    end
end

# Define a custom show method for SIRguided
function Base.show(io::IO, m::SIRguided{T, S}) where {T, S}
    println(io, "SIRguided Process:")
    println(io, "  ξ (1-ξ is probability of new exogeneous infection): ", m.ξ)
    println(io, "  λ (infection rate): ", m.λ)
    println(io, "  μ (recovery rate): ", m.μ)
    println(io, "  ν (rate to get susceptible again): ", m.ν)
    println(io, "  τ (discretisation time step): ", m.τ)
    println(io, "𝒩, ℐ, 𝒪 and O not displayed")
    #println(io, "  𝒩 (neighbors structure): ")
    # for (i, neighbors) in enumerate(m.𝒩)
    #     println(io, "    Neighbor $i: ", neighbors)
    # end
    # println(io, "  ℐ (infected neighbors at each time not displayed) ")
    # for (i, infected_neighbors) in enumerate(m.ℐ)
    #     println(io, "    Time $i: ", infected_neighbors)
    # end
    # println(io, "  𝒪 (some vector): ", m.𝒪)
    # println(io, "  O (3x3 SMatrix): \n", m.O)
end



"""
    initialise_infected_neighbours(𝒪)

    Very simplest way to guess nr of infected neighbours by simply counting the fraction of infections among all observations
    Then for any particles and any time, this nr is assigned 
"""
function initialise_infected_neighbours(𝒪)
    n_times = length(𝒪)
    n_particles= length(𝒪[1].x)
    Xobs = [O.x for O in 𝒪]
    Xobs_flat = vcat(Xobs...)
    frac_infected_observed = sum(Xobs_flat .== _I_)/(length(Xobs_flat) - sum(Xobs_flat .== _L_))
    [fill(frac_infected_observed, n_particles) for _ ∈ 1:n_times] # of course these obs schemes use some bias but fine if only first step
end 


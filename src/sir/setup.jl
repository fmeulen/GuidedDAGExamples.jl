# the following code is taken from setup.jl
# https://github.com/fmeulen/BackwardFilteringSIR/blob/main/setup.jl

# Template transition model constructor
function ptemplate(i, parentstates::Tuple{Vararg{State}}, θ, δ, τ)
    λ, μ, ν = θ



    neighbours = setdiff(1:length(parentstates), i)
    N_i = sum(parentstates[neighbours] .== _I_)

    ψ(u) = exp(-τ*u)

    parentstates[i] == _S_ ? [(1.0-δ)*ψ(λ*N_i),  1.0-(1.0-δ)*ψ(λ*N_i),  0.      ] :
    parentstates[i] == _I_ ? [ 0.,               ψ(μ),                  1.0-ψ(μ)] :
    parentstates[i] == _R_ ? [ 1.0-ψ(ν),         0.,                    ψ(ν)    ] : println("Error")
end


# The general state-space
E = [_S_, _I_, _R_]
statespace = Dict([i => E for i in 1:N])

parents = Dict(i => intersect(i-size_neighbourhood:i+size_neighbourhood, 1:N) for i in 1:N)

nodeToTypestart  = Dict(i => string(i) for i = 1:size_neighbourhood)
nodeToTypemiddle = Dict(c => "Inner" for c ∈ 1+size_neighbourhood:N-size_neighbourhood)
nodeToTypeend    = Dict(i => string(i) for i = N+1-size_neighbourhood:N)

nodeToType = merge(nodeToTypestart, nodeToTypemiddle, nodeToTypeend)

## type To P

function pLeft(parentstates::NTuple{s, State}, θ, δ, τ) where s
    ptemplate(s-size_neighbourhood, parentstates, θ, δ, τ)
end

function pInner(parentstates::NTuple{1+2*size_neighbourhood}, θ, δ, τ)
    ptemplate(1+size_neighbourhood, parentstates, θ, δ, τ)
end

function pRight(parentstates::NTuple{s, State}, θ, δ, τ) where s
    ptemplate(1+size_neighbourhood, parentstates, θ, δ, τ)
end

typeToPstart = Dict(string(i) => pLeft for i = 1:size_neighbourhood)
typeToPinner = Dict("Inner" => pInner)
typeToPend   = Dict(string(i) => pRight for i=N-size_neighbourhood+1:N)

typeToP = merge(typeToPstart, typeToPinner, typeToPend)


## type To Support

typeToSupportstart  = Dict([string(i) => ntuple(j -> E, i+size_neighbourhood) for i=1:size_neighbourhood])
typeToSupportmiddle = Dict(["Inner" => ntuple(j -> E, 1+2*size_neighbourhood)])
typeToSupportend    = Dict([string(i) => ntuple(j -> E, N-i+1+size_neighbourhood) for i = N+1-size_neighbourhood:N])

typeToSupport = merge(typeToSupportstart, typeToSupportmiddle, typeToSupportend)

## done

dynamics(θ, δ, τ) = cpds(nodeToType, typeToP, typeToSupport, θ, δ, τ)

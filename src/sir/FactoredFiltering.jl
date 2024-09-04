# the following code is taken from FactoredFiltering.jl
# https://github.com/fmeulen/BackwardFilteringSIR/blob/main/FactoredFiltering.jl

struct Message{T}
    factoredhtransform::Dict{Int, Vector{T}}
    approximatepullback::Dict{Int, Vector{T}}
end

struct FactorisedMarkovChain{T}
    statespace::Dict{Int, <: Vector{<: Enum}}
    parents::Dict{Int, UnitRange{Int}}
    cpd::Dict{Int, Dict{Vector{UInt8}, Vector{T}}}
    cpdmatrix::Dict{Int, SArray{S, T, 2, L} where {S <: Tuple, L}}
    root::Vector{UInt8}
    N::Int
    T::Int
end

function FactorisedMarkovChain(statespace, parents, (cpd, cpdmatrix), root, (N, T))
    FactorisedMarkovChain(statespace, parents, cpd, cpdmatrix, UInt8.(root), N, T)
end

function inplacemultiplication!(a::AbstractVector, b::AbstractVector)
    for i in eachindex(a)
        @inbounds a[i]=a[i]*b[i]
    end
end

function backwardfiltering(FMC::FactorisedMarkovChain{T}, kernel::Function, approximations, observations, Πroot, size_neighbourhood) where T
    # we first initialise the (faactored) htransforms for all time steps

    htransforms = Dict([t => Dict([i => Vector{T}(ones(length(E))) for i in 1:FMC.N]) for t in 1:FMC.T])
    messages = Dict{Int, Message{T}}()

    obsparents, obscpds, obsstate = observations

    # then we write the observations into the correct htransforms
    # this seems a bit hack-y now but made sense at the time (pullback from individual nodes, obs, etc..)

    for ((i, t), obs) in pairs(obsstate)
        chtransform = obscpds[(i,t)][:,obs]
        inplacemultiplication!(htransforms[t][i], chtransform)
    end

    # Normalise guiding term at time step T

    for i in 1:FMC.N
        htransforms[FMC.T][i] = htransforms[FMC.T][i] / sum(htransforms[FMC.T][i])
    end


    for t in FMC.T:-1:2
        # first computea all the pullback factors. these are passed to a kernel function which produces
        # some projection of the pullback onto the space of marginaals. in this case we use a junction tree algo / boyen koller / fully factorised EP

        factoredpullback = Dict([i => Vector(FMC.cpdmatrix[i]*htransforms[t][i]) for i in 1:FMC.N])

        # this 'approximations' variable was originally the 'reference measures' with which we perform backward marginalisation
        # now it is (ab)used to pass the observation-pullbacks to the junction tree algo -> very hack-y...
        # just always set approxiations to false

        if approximations == false
            approximatepullback = kernel(factoredpullback, htransforms[t-1], FMC.N, t, size_neighbourhood)
        else
            approximatepullback = kernel(factoredpullback, approximations[t-1], FMC.N, t, size_neighbourhood)
        end

        # note: the kernel function handles normalisation and multiplication with observation term

        for i in 1:FMC.N
            #inplacemultiplication!(htransforms[t-1][i], approximatepullback[i])
            htransforms[t-1][i] = approximatepullback[i]
        end

        messages[t] = Message(htransforms[t], approximatepullback)
    end

    logh = 0.0
    for i=1:FMC.N
        logh += log( dot( Πroot[i],  htransforms[1][i] ) )
     end

#     messages[1] = Message(htransforms[1] , htransforms[1])

### what M and F coded up
    # logh = 0.0   # only relevant if θ is updated
    # for (key, value) in htransforms[2]
    #     #a = htransforms[2][key] .* Πroot
    #     a = dot(htransforms[2][key],  Πroot[key])
    #     htransforms[1][key] = [a]
    #     logh += log(a)
    # end

    # messages[1] = Message(htransforms[1] , htransforms[2])
###


#    logh = sum(log(htransforms[1][i][FMC.root[i]]) for i=1:FMC.N) # should not be FMC.root but x0
    messages, logh
end

function cpds(nodeToType, typeToP, typeToSupport, θ::Vector{T}, δ, τ) where T
    typeToCpd = Dict{String, Dict{Vector{UInt8}, Vector{T}}}()
    typeToK = Dict{String, SArray{S, T, 2, L} where {S <: Tuple, L}}()

    nodeToCpd = Dict{Int, Dict{Vector{UInt8}, Vector{T}}}()
    nodeToK = Dict{Int, SArray{S, T, 2, L} where {S <: Tuple, L}}()

    for (type, p) in typeToP
        # should build the dict and matrix in one loop // make this a view
        tempcpd = Dict([UInt8.(collect(parentsstate)) => p(parentsstate, θ, δ, τ) for parentsstate in Iterators.product(typeToSupport[type]...)])
        typeToCpd[type] = tempcpd
        typeToK[type] = SMatrix{length(E)^length(typeToSupport[type]), length(E), T}(reduce(hcat, [tempcpd[UInt8.(collect(parentsstate))] for parentsstate in Iterators.product(typeToSupport[type]...)])')
    end
    # re: memory save: assigning nodetoK keys to typetoK values does NOT copy K, only reference
    for (node, type) in nodeToType                                               # potential memory save
        nodeToCpd[node] = typeToCpd[type]                                        # nodeToCpd(node) = typeToCpd[nodeToType[node]]
        nodeToK[node] = typeToK[type]                                            # nodeToK(node) = typeToK[nodeToType[node]]
    end

    nodeToCpd, nodeToK
end

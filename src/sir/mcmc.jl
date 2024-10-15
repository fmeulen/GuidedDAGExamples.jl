innovations(n_times, n_particles) = [rand(n_particles) for _ in 1:n_times]

"""
    update!(Zᵒ, Z, δ, block)

    update innovations Zᴼ for all times in block, Z remains unchanged
    each element of Zᵒ gets written into
"""
function update!(Zᵒ, Z, δ, block)
    ℒ = Uniform(-δ, δ)
    for i in 1:length(Z[1]), t in 1:length(Z)
        Zᵒ[t][i] = t ∈ block ? mod(Z[t][i] + rand(ℒ), 1) : Z[t][i]
    end
end

function mcmc_iteration!(P, Π, Z, Zᵒ, X, Xᵒ, B, ll, prior, blocks, δ, i, printskip)
    acc = 0
    for block in blocks
        update!(Zᵒ, Z, δ, block)
        llᵒ = forward!(Xᵒ, P, Π, B, Zᵒ, prior)
        if log(rand()) < llᵒ - ll
            mod(i, printskip) == 0 && println(i, "  ", ll, "  ", llᵒ, "  ", llᵒ - ll, "  accepted")
            ll = llᵒ
            Z, Zᵒ = Zᵒ, Z
            X, Xᵒ = Xᵒ, X
            acc += 1
        else
            mod(i, printskip) == 0 && println(i, "   ", ll, "  ", llᵒ, "  ", llᵒ - ll, "  rejected")
        end
    end
    return ll, Z, Zᵒ, X, Xᵒ, acc
end





function updatepar!(X, Xᵒ, P, Pᵒ, B, Bᵒ, Π, Z, prior, ll, logproposalratios)
    backward!(Bᵒ, Pᵒ)
    llᵒ = forward!(Xᵒ, Pᵒ, Π, Bᵒ, Z, prior)

    logA = llᵒ - ll  + logproposalratios
    if log(rand()) < logA
        ll = llᵒ
        P, Pᵒ = Pᵒ, P
        X, Xᵒ = Xᵒ, X
        B, Bᵒ = Bᵒ, B
     end
    ll, P, Pᵒ, X, Xᵒ, B, Bᵒ
end

function update_parameters!(P, X, Xᵒ, B, Bᵒ, Π, Z, ll, propσ, prior)
    # Update μ
#     μᵒ = P.μ * exp(propσ * randn())
#     logproposalratios = log(μᵒ/P.μ)
#     Pᵒ = @set P.μ = μᵒ
#     #ll, P, X, Xᵒ, B = updatepar!(X, Xᵒ, P, Pᵒ, B, Bᵒ, Π, Z, prior, ll, logproposalratios)
#     ll, P, Pᵒ, X, Xᵒ, B, Bᵒ = updatepar!(X, Xᵒ, P, Pᵒ, B, Bᵒ, Π, Z, prior, ll, logproposalratios)

#     # Update λ
#     λᵒ = P.λ * exp(propσ * randn())
#     logproposalratios = log(λᵒ/P.λ)
#     Pᵒ = @set P.λ = λᵒ
# #    ll, P, X, Xᵒ, B = updatepar!(Xᵒ, X, Pᵒ, P, Π, B, Z, prior, ll, logproposalratios)
#     ll, P, Pᵒ, X, Xᵒ, B, Bᵒ = updatepar!(X, Xᵒ, P, Pᵒ, B, Bᵒ, Π, Z, prior, ll, logproposalratios)

    # Update ν
    νᵒ = P.ν * exp(propσ * randn())
    logproposalratios = log(νᵒ/P.ν) 
    Pᵒ = @set P.ν = νᵒ
    #ll, P, X, Xᵒ, B = updatepar!(Xᵒ, X, Pᵒ, P, Π, B, Z, prior, ll, logproposalratios)
    ll, P, Pᵒ, X, Xᵒ, B, Bᵒ = updatepar!(X, Xᵒ, P, Pᵒ, B, Bᵒ, Π, Z, prior, ll, logproposalratios)

    #return ll, P, X, Xᵒ, B
    ll, P, X, Xᵒ, B, Bᵒ
end



# function mcmc(P::SIRguided, Π, Z, blocks;  δ = 0.1, γ = 0.7,
#     acc = 0, ITER = 100, adaptmax=1000,
#     par_estimation = true,
#     propσ = 0.1,
#     prior = (μ=Exponential(5.0), λ = Exponential(5.0), ν=Exponential(5.0)),
#     printskip=5)
    

#    @unpack 𝒪 = P
#    n_times, n_particles = length(𝒪), length(𝒪[1].x)
#    n_blocks = length(blocks)

#    B = backward(P)
#    #Z = innovations(n_times, n_particles)
#    Zᵒ = deepcopy(Z)
#    X, ll  = forward(P, Π, B, Z)
#    Xᵒ = deepcopy(X)

#    XX = [copy(X)]
#    lls = [ll]
#    θs = [param(P)]

#    for i in 1:ITER
#        for block in blocks
#            update!(Zᵒ, Z, δ, block)
#            llᵒ = forward!(Xᵒ, P, Π, B, Zᵒ)

#            if log(rand()) < llᵒ - ll
#                mod(i,printskip)==0 && println(i,"  ",ll,"  ", llᵒ,"  ", llᵒ-ll, "  accepted")
#                ll = llᵒ
#                for t ∈ block
#                    for i in 1:n_particles
#                        Z[t][i] = Zᵒ[t][i] # you can just swap the objects by reference i think
#                    end
#                end

#                for t in 1:n_times
#                    for i in 1:n_particles
#                        X[t][i] = Xᵒ[t][i]
#                    end
#                end
#                acc += 1
#            else
#                mod(i,printskip)==0 && println(i, "   ", ll,"  ", llᵒ,"  ", llᵒ-ll, "  rejected")
#            end
#            i ÷ 10 == 0 && push!(XX, deepcopy(X))
#            push!(lls, ll)
#        end
#        if i < adaptmax && i÷10==0
#            ℐnew = γ * P.ℐ  + (1.0-γ) * count_infections(X, 𝒩)
#            @reset P.ℐ = ℐnew
#            B = backward(P)
#            ll = forward!(X, P, Π, B, Z)
#        end

#        if par_estimation
#            # update μ
#            μᵒ = P.μ * exp(propσ*randn())
#            logprior_proposalratios = (log(μᵒ) - log(P.μ)) + logpdf(prior.μ,μᵒ) - logpdf(prior.μ,P.μ)
#            Pᵒ = @set P.μ = μᵒ
#            ll, P, B = updatepar!(Xᵒ, X, Pᵒ, P, Π,  B, Z,  ll, logprior_proposalratios)

#            # update λ
#            λᵒ = P.λ * exp(propσ*randn())
#            logprior_proposalratios = (log(λᵒ) - log(P.λ)) + logpdf(prior.λ,λᵒ) - logpdf(prior.λ,P.λ)
#            Pᵒ = @set P.λ = λᵒ
#            ll, P, B = updatepar!(Xᵒ, X, Pᵒ, P, Π,  B, Z,  ll, logprior_proposalratios)


#            # update ν
#            νᵒ = P.ν * exp(propσ*randn())
#            logprior_proposalratios = (log(νᵒ) - log(P.ν)) + logpdf(prior.ν,νᵒ) - logpdf(prior.ν,P.ν)
#            Pᵒ = @set P.ν = νᵒ
#            ll, P, B = updatepar!(Xᵒ, X, Pᵒ, P, Π,  B, Z, ll, logprior_proposalratios)
#            push!(θs, param(P))
#        end
#    end
#    @show acc/(ITER*n_blocks)
#    XX, lls, θs, P
# end


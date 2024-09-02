innovations(n_times, n_particles) = [rand(n_particles) for _ in 1:n_times]

"""
    update!(Zᵒ, Z, δ, block)

    update innovations Zᴼ for all times in block
"""
function update!(Zᵒ, Z, δ, block)
    n_particles = length(Z[1])
    n_times = length(Z)
    ℒ = Uniform(-δ, δ)
    for i in 1:n_particles
        for t ∈ 1:n_times         
            Zᵒ[t][i] = t ∈ block ? mod(Z[t][i] + rand(ℒ), 1) : Z[t][i]
        end
    end
end


function mcmc(𝒪, P::SIRguided, Π, blocks;  δ = 0.1, γ = 0.7,
     acc = 0, ITER = 100, adaptmax=1000, 
     par_estimation = true,
     propσ = 0.1, 
     prior = (μ=Exponential(5.0), λ = Exponential(5.0), ν=Exponential(5.0))
     )


    n_times, n_particles = length(𝒪), length(𝒪[1].x)
    n_blocks = length(blocks)
    
    B, ll0 = backward(P, 𝒪)
    Z = innovations(n_times, n_particles)
    Zᵒ = deepcopy(Z)
    X, ll  = forward(P, Π, B, Z, ll0)
    Xᵒ = deepcopy(X)

    XX = [copy(X)]
    lls = [ll]
    θs = [param(P)]

    for i in 1:ITER
        for block in blocks
            update!(Zᵒ, Z, δ, block)
            llᵒ = forward!(Xᵒ, P, Π, B, Zᵒ, ll0)
                          
            if log(rand()) < llᵒ - ll
                mod(i,500)==0 && println(i,"  ",ll,"  ", llᵒ,"  ", llᵒ-ll, "  accepted")
                ll = llᵒ
                for t ∈ block
                    for i in 1:n_particles
                        Z[t][i] = Zᵒ[t][i] 
                    end
                end

                for t in 1:n_times
                    for i in 1:n_particles
                        X[t][i] = Xᵒ[t][i] 
                    end
                end
                acc += 1
            else 
                mod(i,500)==0 && println(i, "   ", ll,"  ", llᵒ,"  ", llᵒ-ll, "  rejected")
            end
            i ÷ 10 == 0 && push!(XX, deepcopy(X))
            push!(lls, ll)
        end
        if i < adaptmax && i÷10==0
            ℐnew = γ * P.ℐ  + (1.0-γ) * count_infections(X, 𝒩)
            @reset P.ℐ = ℐnew
            B, ll0 = backward(P, 𝒪)
            ll = forward!(X, P, Π, B, Z, ll0)
        end
    
        if par_estimation
            # update μ
            μᵒ = P.μ * exp(propσ*randn())
            logprior_proposalratios = (log(μᵒ) - log(P.μ)) + logpdf(prior.μ,μᵒ) - logpdf(prior.μ,P.μ)
            Pᵒ = @set P.μ = μᵒ
            ll, ll0, P, B = updatepar!(Xᵒ, X, Pᵒ, P, Π, 𝒪, B, Z, ll0, ll, logprior_proposalratios)

            # update λ
            λᵒ = P.λ * exp(propσ*randn())
            logprior_proposalratios = (log(λᵒ) - log(P.λ)) + logpdf(prior.λ,λᵒ) - logpdf(prior.λ,P.λ)
            Pᵒ = @set P.λ = λᵒ
            ll, ll0, P, B = updatepar!(Xᵒ, X, Pᵒ, P, Π, 𝒪, B, Z, ll0, ll, logprior_proposalratios)


            # update ν
            νᵒ = P.ν * exp(propσ*randn())
            logprior_proposalratios = (log(νᵒ) - log(P.ν)) + logpdf(prior.ν,νᵒ) - logpdf(prior.ν,P.ν)
            Pᵒ = @set P.ν = νᵒ
            ll, ll0, P, B = updatepar!(Xᵒ, X, Pᵒ, P, Π, 𝒪, B, Z, ll0, ll, logprior_proposalratios)
            push!(θs, param(P))
        end
    end
    @show acc/(ITER*n_blocks)
    XX, lls, θs, P
end




function updatepar!(Xᵒ, X, Pᵒ, P, Π, 𝒪, B, Z, ll0, ll, logprior_proposalratios)
    Bᵒ, ll0ᵒ = backward(Pᵒ, 𝒪)
    llᵒ = forward!(Xᵒ, Pᵒ, Π, Bᵒ, Z, ll0ᵒ)

    logA = llᵒ - ll  + logprior_proposalratios
    if log(rand()) < logA
        ll = llᵒ
        ll0 = ll0ᵒ
        P = Pᵒ
        B, Bᵒ = Bᵒ, B
        for t in eachindex(X)
            for i in eachindex(X[t])
               X[t][i] = Xᵒ[t][i]
            end
        end
    end
    ll, ll0, P, B
end




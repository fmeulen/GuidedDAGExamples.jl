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

# following changes X and B
function updateparλ!(X, P, Π, 𝒪, B, Z, ll)
    propσ = 0.1
    prior = [Exponential(5.0)]
    @unpack λ, μ, ν = P
    λᵒ = λ * exp(propσ*randn())
    Pᵒ = SIRguided(P.ξ, λᵒ, P.μ, P.ν, P.τ, P.𝒩, P.ℐ)
    Bᵒ, logwᵒ = backward(Pᵒ, 𝒪)
    Xᵒ, llᵒ = forward(Pᵒ, Π, Bᵒ, Z, logwᵒ)
    if log(rand()) < llᵒ - ll  + (log(λᵒ) - log(λ)) + logpdf(prior[1],λᵒ) - logpdf(prior[1],λ)
#        if mod(it, skip_print)==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    λ ", λᵒ)  end
        ll, llᵒ, P, Pᵒ = llᵒ, ll, Pᵒ, P
        for t in eachindex(X)
            for i in eachindex(X[1])
               X[t][i] = Xᵒ[t][i]
               B[t][i] = Bᵒ[t][i]
            end
        end
#        accpar[1] += 1
    else
#        if mod(it, skip_print )==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    λ ", P.λ)  end
    end
    ll, P
end

function updateparμ!(Xᵒ, X, P, Π, 𝒪, B, Z, ll)
    propσ = 0.1
    prior = [Exponential(5.0)]
    @unpack λ, μ, ν = P

    μᵒ = μ * exp(propσ*randn())
    Pᵒ = SIRguided(P.ξ, P.λ, μᵒ, P.ν, P.τ, P.𝒩, P.ℐ)
    Bᵒ, logwᵒ = backward(Pᵒ, 𝒪)
    llᵒ = forward!(Xᵒ, Pᵒ, Π, Bᵒ, Z, logwᵒ)
    
    if log(rand()) < llᵒ - ll  + (log(μᵒ) - log(μ)) + logpdf(prior[1],μᵒ) - logpdf(prior[1],μ)
#        if mod(it, skip_print)==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    λ ", λᵒ)  end
#         ll, llᵒ, P, P = llᵒ, ll, Pᵒ, P
#         for t in eachindex(X)
#             for i in eachindex(X[1])
#                X[t][i] = Xᵒ[t][i]
#                B[t][i] = Bᵒ[t][i]
#             end
#         end
# #        accpar[1] += 1
#     else
# #        if mod(it, skip_print )==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    λ ", P.λ)  end
#     end
#     ll, P, B
 #       return(Pᵒ, Bᵒ, llᵒ)
        ll = llᵒ
        Pout = Pᵒ
        #Bout = Bᵒ
        for t in eachindex(X)
            for i in eachindex(X[t])
               X[t][i] = Xᵒ[t][i]
               B[t][i] = Bᵒ[t][i]
            end
        end
    else
  #      return(P, B, ll)
        Pout = P
        #Bout = B
    end
     ll, Pout
end



function mcmc(𝒪, P::SIRguided, Π, blocks;  δ = 0.1, γ = 0.7, acc = 0, ITER = 100)
    adaptmax = 1#ITER÷2

    n_times, n_particles = length(𝒪), length(𝒪[1].x)
    

    B, logw = backward(P, 𝒪)

    Z = innovations(n_times, n_particles)
    Zᵒ = deepcopy(Z)
    X, ll  = forward(P, Π, B, Zᵒ, logw)
    Xᵒ = deepcopy(X)

    XX = [copy(X)]
    lls = [ll]
    θs = [param(P)]

    for i in 1:ITER
        for block in blocks
            update!(Zᵒ, Z, δ, block)
            llᵒ = forward!(Xᵒ, P, Π, B, Zᵒ, logw)
                          
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
        if i < adaptmax
            ℐnew = γ * P.ℐ  + (1.0-γ) * count_infections(X, 𝒩)
            P = SIRguided(P.ξ, P.λ, P.μ, P.ν, P.τ, P.𝒩, ℐnew) 
            # infected_neighbours_new = count_infections(X, 𝒩)
            # infected_neighbours .= γ * infected_neighbours + (1-γ) * infected_neighbours_new
            B, logw = backward(P, 𝒪)
             X, ll  = forward(P, Π, B, Z, logw)
        end

        # par updating
        propσ = 0.1
        prior = [Exponential(5.0)]
        @unpack λ, μ, ν = P
    
        μᵒ = μ * exp(propσ*randn())
        Pᵒ = SIRguided(P.ξ, P.λ, μᵒ, P.ν, P.τ, P.𝒩, P.ℐ)
        Bᵒ, logwᵒ = backward(Pᵒ, 𝒪)
        llᵒ = forward!(Xᵒ, Pᵒ, Π, Bᵒ, Z, logwᵒ)
        
        if log(rand()) < llᵒ - ll  + (log(μᵒ) - log(μ)) + logpdf(prior[1],μᵒ) - logpdf(prior[1],μ)
            ll = llᵒ
            logw = logwᵒ
            μ = μᵒ
            P = SIRguided(P.ξ, P.λ, μ, P.ν, P.τ, P.𝒩, P.ℐ)
            B, Bᵒ = Bᵒ, B
            for t in eachindex(X)
                for i in eachindex(X[t])
                   X[t][i] = Xᵒ[t][i]
                   #B[t][i] = Bᵒ[t][i]
                end
            end
        end


        push!(θs, param(P))

    end
    @show acc/(ITER*n_blocks)
    XX, lls, θs
end

P = SIRguided(Ptrue.ξ, Ptrue.λ, 5.2, Ptrue.ν, Ptrue.τ, Ptrue.𝒩, ℐ)
Xs, lls, θs = mcmc(𝒪, P, Π, blocks; δ=0.1, ITER=10_000);
lo = @layout [a;b]
μs = getindex.(θs,2);
plot(plot(lls), plot(μs), layout=lo)

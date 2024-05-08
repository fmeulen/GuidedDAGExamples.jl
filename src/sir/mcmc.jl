innovations(n_times, n_particles) = [rand(n_particles) for _ in 1:n_times]

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


function mcmc(𝒪, P::SIRguided, Π;  δ = 0.1, γ = 0.7, acc = 0, n_blocks = 4, ITER = 100)
    adaptmax = ITER÷8

    n_times, n_particles = length(𝒪), length(𝒪[1].x)
    blocks = partition_into_blocks_close(n_times, n_blocks)

    # Xobs = [O.x for O in 𝒪]
    # infected_neighbours = count_infections(Xobs, P.𝒩)
    B = backward(P, 𝒪)
    ℱ = forward(P, Π, B)

    Z = innovations(n_times, n_particles)
    Zᵒ = deepcopy(Z)
    X, ll  = ℱ(Z)

    XX = [copy(X)]
    lls = [ll]
    for i in 1:ITER
        for block in blocks
            update!(Zᵒ, Z, δ, block)
            Xᵒ, llᵒ  = ℱ(Zᵒ)
              
            if log(rand()) < llᵒ - ll
                i÷10==0 && println(ll,"  ", llᵒ,"  ", llᵒ-ll, "  accepted")
                ll = llᵒ

                #Z .= Zᵒ
                for t ∈ 1:n_times
                    for i in 1:n_particles
                        Z[t][i] = Zᵒ[t][i] 
                    end
                end

                X .= Xᵒ
                
                acc += 1
            else 
                i÷10==0 && println(ll,"  ", llᵒ,"  ", llᵒ-ll, "  rejected")
            end
            push!(XX, deepcopy(X))
            push!(lls, ll)
        end
        if i < 1#adaptmax
            infected_neighbours_new = count_infections(X, 𝒩)
            infected_neighbours .= γ * infected_neighbours + (1-γ) * infected_neighbours_new
            B = backward(P, 𝒪)
        # X, ll  = forward(P, Π, B, Z)
        end
    end
    @show acc/(ITER*n_blocks)
    XX, lls
end

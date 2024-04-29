innovations(n_times, n_particles) = [rand(n_particles) for _ in 1:n_times]

function update(Z, δ, block)
    n_particles = length(Z[1])
    ℒ = Uniform(-δ, δ)
    Zᵒ = copy(Z)
    Zᵒ[block] = [mod.(Z[t] + rand(ℒ, n_particles), 1) for t in block  ]
    Zᵒ
end


function mcmc(𝒪, P::SIRguided, Π;  δ = 0.1, γ = 0.7, acc = 0, n_blocks = 4, ITER = 100)
    adaptmax = ITER÷8

    n_times, n_particles = length(𝒪), length(𝒪[1].x)
    blocks = partition_into_blocks_close(n_times, n_blocks)

    # Xobs = [O.x for O in 𝒪]
    # infected_neighbours = count_infections(Xobs, P.𝒩)

    Z = innovations(n_times, n_particles)
    B = backward(P, 𝒪)
    X, ll  = forward(P, Π, B, Z)

    Xs = [X]
    lls = [ll]
    for i in 1:ITER
        for block in blocks
            Zᵒ = update(Z, δ, block)
            Xᵒ, llᵒ  = forward(P, Π, B, Zᵒ)
            @show ll, llᵒ
            if log(rand()) < llᵒ - ll
                ll = llᵒ
                for k in eachindex(Z)
                    Z[k] .= Zᵒ[k]
                end
#                Zᵒ, Z = Z, Zᵒ
                Xᵒ, X = X, Xᵒ
                #Zᵒ .= Z
                #Xᵒ .= X
            #    println("acc")
                acc += 1
            end
            push!(Xs, copy(X))
            push!(lls, copy(ll))
        end
        if i < 1#adaptmax
            infected_neighbours_new = count_infections(X, 𝒩)
            infected_neighbours .= γ * infected_neighbours + (1-γ) * infected_neighbours_new
            B = backward(P, 𝒪)
        # X, ll  = forward(P, Π, B, Z)
        end
    end
    @show acc/(ITER*n_blocks)
    Xs, lls
end

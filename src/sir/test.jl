

function f(x)
    push!(x, 4.0)
    y = sin.(x)
    y
end

x = [1.0, 2.0]
f(x)
x

using ConcreteStructs

@concrete struct  Mystruct
    a
end

UU= Mystruct(Z)



 B, logw = backward(P, 𝒪)

 Π = [SA_F64[0.9, 0.1, 0.0] for _ in 1:n_particles]
 Z = innovations(n_times, n_particles)

 XX, ll = forward(P, Π, B, Z, logw)
 XX2, ll2 = forward(P, Π, B, Z, logw) 

 @btime backward(P, 𝒪);
@code_warntype backward(P, 𝒪)
 @btime forward(P, Π, B, Z, logw);
@code_warntype forward(P, Π, B, Z, logw);
# # deepcopy: 59.959 μs (1134 allocations: 85.78 KiB)

# forward!(XX, P, Π, B, Z, logw);
# forward!(XX, P, Π, B, Z, logw);

# function forward2!(Xs, P::SIRguided, Π, B, Z, logw)
#     n_steps, n_particles = length(Z), length(Π)

#     # sample initial state
#     X = Vector{State}(undef, n_particles)
#     z = Z[1]
#     ll = logw
#     for i in 1:n_particles
#         p = Π[i] .* B[1][i]
#         X[i] = rand𝒳(z[i], p/sum(p))
#         ll += log(sum(p))
#     end
    

#     Xs[1] .= X
#     for t in 2:n_steps
#         lw = guide!(X, P, B[t], Z[t], P.ℐ[t-1])
#         ll += lw
#         Xs[t] .= X
#     end
#     Xs, ll
# end

# t=2
# @btime guide!(XX[1], P, B[t], Z[t], P.ℐ[t-1])
wd = @__DIR__
cd(wd)


using Distributions
PLOT = true
using LinearAlgebra
using DataFrames
using Random
using StaticArrays
#using Revise
using ConcreteStructs
if PLOT
    using RCall
    using Plots
else
    macro rput(args...)
    end
    macro R_str(args...)
    end
    macro layout(args...)
    end

end

include("funcdefs.jl")

figdir = mkpath(joinpath(wd,"figs"))

############## generate data
Random.seed!(3)

n_particles = 20
n_times = 50
𝒩 = set_neighbours(n_particles)
ξ, λ, μ, ν, τ =  1.0, 1.5, 2.0, 5.1, 0.1
Ptrue = SIRforward(ξ, λ, μ, ν, τ, 𝒩)
#x0 = vcat(_I_, fill(_S_,n-2),_I_)
x0 = vcat(fill(_I_,7), fill(_S_,n_particles-7))

Xtrue = sample_trajectory(Ptrue::SIRforward, n_times, x0)
# visualise
PLOT && Plots.heatmap(obs2matrix(Xtrue)' ,xlabel="time", ylabel="individual", 
        title="Forward simulated")

δobs = 0.001
O = SA[1.0-δobs δobs/2.0 δobs/2.0; δobs/2.0 1.0- δobs δobs/2.0; δobs/2.0 δobs/2.0 1-δobs]


@concrete struct Observation
    ind         # at each time, indices of observed particles
    h           # corresponding h-vectors
    x           # values at all times (latent = _L_)
end

function create_data(samplesize, n_times, n_particles, O)
    nrobs_at_each_time = rand(Multinomial(samplesize, n_times))  
    ind_obs = Vector{Int64}[]
    for i in 1:n_times
        k = nrobs_at_each_time[i]
        ids = rand(DiscreteUniform(1,n_particles), k)
        push!(ind_obs, ids)
    end

    Xobs = [fill(_L_, n_particles) for _ in 1:n_times   ]
    h = Vector{SVector{3, Float64}}[]
    for t in 1:n_times
        ht = SVector{3, Float64}[]
        for i in ind_obs[t]
            x = Xtrue[t][i]
            println(x)
            Xobs[t][i] = x
            push!(ht, O * observationmessage(x))
        end
        push!(h, ht)
    end

    [Observation(ind_obs[t], h[t], Xobs[t]) for t in 1:n_times]
end

samplesize = 30
𝒪 = create_data(samplesize, n_times, n_particles, O)

Xobs = [𝒪[i].x for i in eachindex(𝒪)]

# visualise

# construct observation ColorPalette
defaultpalette = palette(cgrad(:default, categorical=true), 3)
# white = RGBA{Float64}(255, 255, 255)
# white = RGBA{Float64}(16, 59, 223, 0.12)
white = RGBA(52, 162, 231, 0.23)

observationcolors = vec(hcat(white, defaultpalette.colors.colors...))
observationpalette = ColorPalette(typeof(defaultpalette.colors)(observationcolors, "", ""))


pobs = heatmap(obs2matrix(Xobs)', xlabel="time", ylabel="individual", 
colorbar=true, color=observationpalette, dps=600, title="observed", background_color_subplot=white)









P = SIRguided(.3, 2.0, 0.8, Ptrue.τ, Ptrue.𝒩) # initialisation

ITER = 10_000
BI = div(ITER,2)
#prior =(Uniform(0,20.0),Uniform(0,20.0),Uniform(0,20.0))
prior = (Exponential(5.0), Exponential(5.0), Exponential(1.0))
@time θθ, Xfinal, N, Q, Xinit, Xmid, acc, accpar, difflr =
    sir_inference(Xobs, P, J; ρ=0.99, prior=prior,
                        ITER=ITER, γ =0.7,  propσ = 0.05)


println()
println("True value and posterior mean: ")
println("λ: ", Ptrue.λ,"     ", round(mean(map(x->x.λ, θθ[BI:end]));digits=3))
println("μ: ", Ptrue.μ,"     ", round(mean(map(x->x.μ, θθ[BI:end]));digits=3))
println("ν: ", Ptrue.ν,"     ", round(mean(map(x->x.ν, θθ[BI:end]));digits=3))
println()
println("Mean acceptance probability for segments: $acc")
println("Mean acceptance probability for parameters: $accpar")

include("rplotting.jl")
#include("juliaplotting.jl")

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
using StatsBase

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

include("createdata.jl")
include("funcdefs.jl")
include("backward.jl")

include("guide.jl")
include("mcmc.jl")
include("partition.jl")

figdir = mkpath(joinpath(wd,"figs"))

############## generate data
Random.seed!(3)

n_particles = 12
n_times = 10
𝒩 = set_neighbours(n_particles, 1)
ξ, λ, μ, ν, τ =  1.0, 1.5, 2.0, 3.1, 0.1
Ptrue = SIRforward(ξ, λ, μ, ν, τ, 𝒩)
#x0 = vcat(_I_, fill(_S_,n-2),_I_)
x0 = vcat(fill(_S_,3), [_I_], fill(_S_,7), [_I_], fill(_S_,n_particles-12))

Xtrue = sample_trajectory(Ptrue::SIRforward, n_times, x0)
# visualise
pforward = Plots.heatmap(obs2matrix(Xtrue)' ,xlabel="time", ylabel="individual", title="Forward simulated")
PLOT && pforward

δobs = 0.001
O = SA[1.0-δobs δobs/2.0 δobs/2.0; δobs/2.0 1.0- δobs δobs/2.0; δobs/2.0 δobs/2.0 1-δobs]

function plotpath(X; name="path") 
    n_particles = length(X[1])
    Xc = copy(X)
    push!(Xc, vcat([_S_, _I_, _R_], fill(_L_, n_particles-3)))
    # construct observation ColorPalette
    defaultpalette = palette(cgrad(:default, categorical=true), 3)
    # white = RGBA{Float64}(255, 255, 255)
    # white = RGBA{Float64}(16, 59, 223, 0.12)
    white = RGBA(52, 162, 231, 0.23)
    observationcolors = vec(hcat(white, defaultpalette.colors.colors...))
    observationpalette = ColorPalette(typeof(defaultpalette.colors)(observationcolors, "", ""))
    p = heatmap(obs2matrix(Xc)', xlabel="time", ylabel="individual", 
    colorbar=true, color=observationpalette, dps=600, title=name, background_color_subplot=white)
    return p
end
samplesize = 20
𝒪 = create_data(samplesize, n_times, n_particles, O)

Xobs = [O.x for O in 𝒪]

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

lo = @layout [a;b]
plot(pforward, pobs, layout=lo)

###############################################################

#P = SIRguided(1.0,.3, 2.0, 0.8, Ptrue.τ, Ptrue.𝒩) # initialisation

# construct guided process from 𝒪 and 𝒩
Xobs = [O.x for O in 𝒪]
ℐ = count_infections(Xobs, 𝒩)
P = SIRguided(Ptrue.ξ, Ptrue.λ,  Ptrue.μ, Ptrue.ν, Ptrue.τ, Ptrue.𝒩, ℐ)

#exp_neighb(P,ave_ninf) = (λ=ave_ninf*P.λ, μ=P.μ, ν=P.ν)


B = backward(P, 𝒪)
    
Π = [[0.99, 0.02, 0.0] for _ in 1:n_particles]


Z = innovations(n_times, n_particles)
X, ll  = forward(P, Π, B, Z)
@show ll

lo = @layout [a;b;c]
ptrue = plotpath(Xtrue; name="true")
pobs = plotpath(Xobs; name="observed")
pguided = plotpath(X; name="guided")
plot(ptrue, pobs, pguided, layout=lo)


Xs, lls = mcmc(𝒪, P, Π; δ=0.1, ITER=1000)
plot(lls)

# anim = @animate for x in Xs
#     plotpath(x)
# end
pforward

plotpath(Xs[end])


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

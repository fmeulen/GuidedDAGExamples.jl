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
using UnPack
using Accessors
using BenchmarkTools

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
include("forward.jl")
include("mcmc.jl")
include("partition.jl")
include("plotting.jl")

figdir = mkpath(joinpath(wd,"figs"))

############## generate data
#Random.seed!(30)

n_particles = 30
n_times = 100
samplesize = (n_times * n_particles)÷20

# set neighbourhood structure
𝒩 = set_neighbours(n_particles, 2)

# set true pars
ξ, λ, μ, ν, τ =  1.0, 3.5, 2.0, 3.1, 0.1
Ptrue = SIRforward(ξ, λ, μ, ν, τ, 𝒩)

# set initial state
x0 = vcat(_I_, fill(_S_,n_particles-2),_I_)
#x0 = [_I_, _S_, _S_, _S_, _S_]
#x0 =  vcat(fill(_S_,3), [_I_], fill(_S_,7), [_I_], fill(_S_,n_particles-12))

Xtrue = sample_trajectory(Ptrue::SIRforward, n_times, x0)
# visualise
pforward = plotpath(Xtrue; name="forward")
PLOT && pforward

# set observation scheme
δobs = 0.001
O = SA[1.0-δobs δobs/2.0 δobs/2.0; δobs/2.0 1.0- δobs δobs/2.0; δobs/2.0 δobs/2.0 1-δobs]


𝒪 = create_data(samplesize, n_times, n_particles, O)



# visualise
Xobs = [O.x for O in 𝒪]
pobs = plotpath(Xobs)

lo = @layout [a;b]
plot(pforward, pobs, layout=lo)

###############################################################


# construct guided process from 𝒪 and 𝒩
Xobs = [O.x for O in 𝒪]
ℐ = count_infections(Xobs, 𝒩)
P = SIRguided(Ptrue.ξ, Ptrue.λ,  Ptrue.μ, Ptrue.ν, Ptrue.τ, Ptrue.𝒩, ℐ)

#exp_neighb(P,ave_ninf) = (λ=ave_ninf*P.λ, μ=P.μ, ν=P.ν)
Xobs_flat = vcat(Xobs...)

# set guided process
frac_infected_observed = sum(Xobs_flat .== _I_)/(length(Xobs_flat) - sum(Xobs_flat .== _L_))
ℐ = [fill(frac_infected_observed, n_particles) for _ ∈ 1:n_times]
P = SIRguided(Ptrue.ξ, Ptrue.λ ,  Ptrue.μ, Ptrue.ν, Ptrue.τ, Ptrue.𝒩, ℐ)

B, logw = backward(P, 𝒪)

# set prior
Π = [SA_F64[0.9, 0.1, 0.0] for _ in 1:n_particles]


############################################################
##### this can go later #######
 Z = innovations(n_times, n_particles)
 X, ll  = forward(P, Π, B, Z, logw);
 ll
loglikelihood(X, Π, B, 𝒪, O) # should be the same as ll
# @show ll

# Y = copy(X);
# ll = forward!(Y, P, Π, B, Z, logw);
# @show ll

# t=2
# @btime guide!(Y[t], P, B[t], Z[t], P.ℐ[t-1])

# lo = @layout [a;b;c]
# ptrue = plotpath(Xtrue; name="true")
# pobs = plotpath(Xobs; name="observed")
# pguided = plotpath(X; name="guided")
# println(ll)
# plot(ptrue, pobs, pguided, layout=lo)
# @show ll

# Zᵒ = deepcopy(Z)
# update!(Zᵒ,    Z, 0.3, 1:4);
# Xᵒ, llᵒ  = forward(P, Π, B, Zᵒ, logw);
# @show ll, llᵒ, llᵒ-ll
############################################################


P = SIRguided(Ptrue.ξ, Ptrue.λ,  Ptrue.μ, Ptrue.ν, Ptrue.τ, Ptrue.𝒩, ℐ)
#P = SIRguided(Ptrue.ξ, 5.0, .2, 7.0, Ptrue.τ, Ptrue.𝒩, ℐ)

blocksize = 4
n_blocks= n_times ÷ blocksize
blocks = make_partition(n_times, n_blocks)
blocks = [1:1, 2:5, 6:n_particles]


Xs, lls, θs, P = mcmc(𝒪, P, Π, blocks; δ=0.1, 
                ITER=20_000, 
                adaptmax= 1000,
                par_estimation=true);#false);
lo = @layout [a;b;c;d]
λs = getindex.(θs,:λ);
μs = getindex.(θs,:μ);
νs = getindex.(θs,:ν);
plot(plot(lls), plot(λs, title="λ", label=""), 
                plot(μs, title="μ", label=""), 
                plot(νs, title="ν", label=""), 
                layout=lo, size=(400, 600))

mean(λs) 
mean(μs) 
mean(νs)
params(Ptrue)

lo = @layout  [a;b;c]
plot(plotpath(Xtrue;name="true"),
plotpath(Xs[end]; name="guided"),
plotpath(Xobs;name="observed"),
layout=lo, size=(700,700))

L = length(Xs)
 anim = @animate for i in eachindex(Xs)
     plot(plotpath(Xtrue;name="true"),
          plotpath(Xs[i]; name="$i of $L"),
          plotpath(Xobs;name="observed"),
          layout=lo, size=(400, 600))
 end

 gif(anim, "anim.gif", fps=5)



## ggplot
using DataFrames
using RCall
len = length(μs)
dout = DataFrame(value=vcat(μs, λs, νs), 
                iteration=repeat(1:len, outer=3), 
                parameter=repeat(["mu", "lambda", "nu"], inner=len))

    @rput dout
    R"""
    library(tidyverse)
    str(dout)
    dout %>% ggplot(aes(x=iteration, y=value)) + 
            geom_path() + 
            facet_wrap(~parameter, ncol=1)
    """





prior = (Exponential(5.0), Exponential(5.0), Exponential(1.0))


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

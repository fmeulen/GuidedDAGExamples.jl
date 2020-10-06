using Distributions
PLOT = true
using LinearAlgebra
using DataFrames
using Random
using StaticArrays
using Revise
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

wd = @__DIR__
cd(wd)
include("funcdefs.jl")

figdir = mkpath(joinpath(wd,"figs"))

############## generate data
Random.seed!(3)
#Random.seed!(9)
n = 100
𝒩 = set_neighbours(n)
λ= 2.5; μ=0.6; ν=0.1; τ = 0.1
Ptrue = SIRforward(λ, μ, ν, τ, 𝒩)
#x0 = vcat(_I_, fill(_S_,n-2),_I_)
x0 = vcat(fill(_I_,7), fill(_S_,n-7))


# generate discrete time observations
J = 50 # then on each segment we impute J-2 points, because sample segment returns J values, including start and endpoint
nobs = 10
tobs = 1:(J-1):nobs*(J-1)
collect(tobs)
tfull = 1:tobs[end]
Xtrue = sample(Ptrue::SIRforward,tfull[end],x0)
Xobs = Xtrue[tobs]

# visualise
PLOT && Plots.heatmap(1:n, tfull, obs2matrix(Xtrue),title="Forward simulated")
PLOT && Plots.heatmap(1:n, tobs, obs2matrix(Xobs),title="Observed")

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

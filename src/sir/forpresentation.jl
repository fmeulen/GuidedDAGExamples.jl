wd = @__DIR__
cd(wd)

figdir = mkpath(joinpath(wd,"figs"))

using Distributions
#PLOT = true
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

using Plots


include("funcdefs.jl")
include("createdata.jl")
include("backward.jl")
include("forward.jl")
include("mcmc.jl")
include("partition.jl")
include("plotting.jl")
include("filtering.jl")


## make simple plots to explain the dynamics, no random numbers
τ = 0.1 # time step

prior = (μ=Exponential(5.0), λ = Exponential(5.0), ν=Exponential(5.0))

presfigdir = mkpath(joinpath(wd,"figs","presentation"))

n_particles = 6
n_times = 7

latent = fill(fill(_L_, n_particles),n_times)

X =  [
[_S_, _S_, _I_, _S_,_S_,_S_],
[_S_, _S_, _I_, _S_,_S_,_S_],
[_S_, _I_, _I_, _S_,_S_,_S_],
[_I_, _I_, _R_, _S_,_S_,_S_],
[_I_, _I_, _R_, _S_,_S_,_S_],
[_I_, _I_, _S_, _S_,_S_,_S_],
[_I_, _R_, _S_, _S_,_S_,_S_]
]

plotpath(X)
for i in 1:n_times
    Xpart = copy(latent)
    Xpart[1:i] = X[1:i]
    p = plotpath(Xpart;name="")
    png(p, presfigdir*"/smallexample$i.png")
end

anim = @animate for i in 1:n_times
    Xpart = latent
    Xpart[1:i] = X[1:i]
    p = plotpath(Xpart;name="")
    p
end
mp4(anim,presfigdir*"/smallexample.mp4", fps=1.5)




## now a bigger example 
n_particles = 50
n_times = 100
size_neighbourhood = 2
latent = fill(fill(_L_, n_particles),n_times)

Random.seed!(30)

# set neighbourhood structure
𝒩 = set_neighbours(n_particles, size_neighbourhood)

# set true pars
ξ, λ, μ, ν = 0.999, 4.0, 2.5, 0.5
Ptrue1 = SIRforward(ξ, λ, μ, ν, τ, 𝒩)

ξ, λ, μ, ν = 1.0, 4.0, 2.5, 1.5
Ptrue2 = SIRforward(ξ, λ, μ, ν, τ, 𝒩)

# set initial state
x0 = vcat(_I_, fill(_S_,n_particles-2),_I_) 

X = sample_trajectory(Ptrue1, n_times, x0)
plotpath(X)
anim2a = @animate for i in 1:n_times
    Xpart = latent
    Xpart[1:i] = X[1:i]
    p = plotpath(Xpart)
    p
end
mp4(anim2a,presfigdir*"/large_example.mp4", fps=10)

# What if we would run it for way longer?
#Ptrue = SIRforward(ξ, 2.0, 1.0, 0.1, τ, 𝒩) # reduce μ, so stay longer infected, don't get infected easily, reduce λ
n_times_long = 1000
latent = fill(fill(_L_, n_particles),n_times_long)

Random.seed!(30)
X = sample_trajectory(Ptrue1, n_times_long, x0)
plotpath(X)
anim2a_long = @animate for i in 1:n_times_long
    Xpart = latent
    Xpart[1:i] = X[1:i]
    p = plotpath(Xpart)
    p
end
mp4(anim2a_long,presfigdir*"/large_example_long.mp4", fps=200)



# set true pars back to initial settings

# other realisations
Random.seed!(2)
for i in 1:8
    Xtrue = sample_trajectory(Ptrue1, n_times, x0)
    pforward = plotpath(Xtrue; name="")
    png(pforward, presfigdir*"/large_example$i.png")
end 



# also when higher rate of getting again susceptbile

Random.seed!(2)
for i in 1:8
    Xtrue = sample_trajectory(Ptrue2, n_times, x0)
    pforward = plotpath(Xtrue; name="")
    png(pforward, presfigdir*"/large_example_faster_susceptible$i.png")
end 


#################### Generating partical observations #################################

# set observation scheme
case = 2 
if case==1
    Ptrue = Ptrue1
    δobs = 0.0 
    Random.seed!(666) # nice!!
end
if case == 2   
    Ptrue = Ptrue2
    δobs = 0.001
    Random.seed!(79)
end

O = SA[1.0-δobs δobs/2.0 δobs/2.0; δobs/2.0 1.0- δobs δobs/2.0; δobs/2.0 δobs/2.0 1-δobs]

# take some more particles 
n_particles = 100
𝒩large = set_neighbours(n_particles, size_neighbourhood)
Ptrue = @set Ptrue.𝒩 = 𝒩large

n_times = 300

x0 = vcat(_I_, fill(_S_,n_particles-2),_I_) 

# simulate forward path 
Xtrue = sample_trajectory(Ptrue, n_times, x0)  
plotpath(Xtrue)


samplesize = (n_times * n_particles)÷25
𝒪1 = create_data(Xtrue, samplesize, n_times, n_particles, O)

obs_times = 25:25:n_times
𝒪2 = create_data_regular(Xtrue, obs_times, n_times, n_particles, O)

obs_particles = 5:7:n_particles 
𝒪3 = create_data_see_some_particles(Xtrue, obs_particles, n_times, n_particles, O)

#



# visualise
#lo = @layout [a;b]
pforward = plotpath(Xtrue;name="True, unobserved")
png(pforward,  presfigdir*"/large_example_forward_$case.png")

#pforward.attr[:size]
pobs = plotpath(𝒪1;name="What we observe")
png(pobs,  presfigdir*"/large_example_observe1_$case.png")
#pp = plot(pforward, pobs, layout=lo, margi)
#png(pp,  presfigdir*"/large_example_forw_and_observe1.png")

pobs2 = plotpath(𝒪2;name="What we observe")
png(pobs2,  presfigdir*"/large_example_observe2_$case.png")
# pp = plot(pforward, pobs2, layout=lo, size=(600,800))
# png(pp,  presfigdir*"/large_example_forw_and_observe2.png")

pobs3 = plotpath(𝒪3;name="What we observe")
png(pobs,  presfigdir*"/large_example_observe3_$case.png")
# pp = plot(pforward, pobs3, layout=lo)
# png(pp,  presfigdir*"/large_example_forw_and_observe3.png")


#################### GUIDING #################################
Random.seed!(58) 


# choose observations scheme
𝒪 = 𝒪2

# Omessages determines the message sent from observations
δ = δobs # in the guided process assume some noise 
Omessages = SA[1.0-δ δ/2.0 δ/2.0; δ/2.0 1.0- δ δ/2.0; δ/2.0 δ/2.0 1-δ]

pobs = plotpath(𝒪;name="What we observe")
ℐ =  initialise_infected_neighbours(𝒪)
P = SIRguided(Ptrue, ℐ, 𝒪, Omessages) 
Πᵢ = 2.0/n_particles; Πₛ = 1.0 - Πᵢ
Π = [SA_F64[Πₛ, Πᵢ, 0.0] for _ in 1:n_particles]
B = backward(P);


# sample guided process 
Z = innovations(n_times, n_particles)
@time Xguided, ll  = forward(P, Π, B, Z, prior);
ll_rounded = round(ll; digits=1)
pguided = plotpath(Xguided; name="Reconstructed, ll is $ll_rounded")

lo = @layout [a;b;c]
pp = plot(pforward, pguided, pobs, layout=lo)
png(pp,  presfigdir*"/large_example_forward_guided_$case.png")

# sample multiple guided processes

# Zbest = innovations(n_times, n_particles)
# Xguidedbest, llbest = forward(P, Π, B, Zbest, prior)

# lo = @layout [a;b]
# anim4 = @animate for i in 1:30
# global llbest    
#     Z = innovations(n_times, n_particles)
#     Xguided, ll  = forward(P, Π, B, Z, prior);
#     lll =  round(ll; digits=1)
#     pguided = plotpath(Xguided; name="Reconstructed. $lll")
#     if ll > llbest
#         Zbest .= Z
#         llbest = ll
#     end
#     @show ll
#     pp = plot(pforward, pguided, layout=lo)#, size=(400, 790))
#     pp
# end
# @show llbest

# mp4(anim4,presfigdir*"/multipleguided.mp4", fps=2)


# Xguidedbest, llbest = forward(P, Π, B, Zbest, prior)
# lo = @layout [a;b;c]
# plot(pforward, plotpath(Xguidedbest), pobs, layout=lo)


# perhaps do smc lagged filtering from here?

n_ensemble = 10_000
# filter_times = [0, 38, 40] # should include 0
filter_times =  vcat(0,100:100:ntimes(P))
#filter_times = vcat(0, obstimes(P))

Random.seed!(2)
𝕏, logweights = filter(P, Π, prior, n_ensemble, filter_times);

plotpath(𝕏[1]; xlims_=(1,n_times))
plotpath(𝕏[3]; xlims_=(1,n_times))
plotpath(𝕏[4]; xlims_=(1,n_times))

p = plotpath(𝕏[1]; name="Reconstructed")
png(p,  presfigdir*"/large_example_reconstructed2_$case.png")

# so we have learned roughly where the infections are 
ℐ_up =count_infections(𝕏[1], P.𝒩)
P_up = SIRguided(Ptrue, ℐ_up, 𝒪, Omessages) 

𝕏_up, logweights_up = filter(P_up, Π, prior, n_ensemble, filter_times);
p_up = plotpath(𝕏_up[1]; name="Reconstructed")
png(p_up,  presfigdir*"/large_example_reconstructed2up_$case.png")


plot(p, p_up)

# alternatively backward filter all (this can be bad)
Z = innovations(ntimes(P), nparticles(P))
Xguided, ll  = forward(P_up, Π, B, Z, prior);
ll_rounded = round(ll; digits=1)
pguided = plotpath(Xguided; name="Reconstructed, ll is $ll_rounded")

lo = @layout [a;b;c]
pp = plot(pforward, pguided, pobs, layout=lo)


# anim_smc = @animate for  x ∈ 𝕏
#     #plotpath(𝕏[1]; xlims_=(1,n_times))
#     plot(pforward, plotpath(x), pobs, layout=lo)
# end

# mp4(anim_smc,presfigdir*"/smc_example_$case.mp4", fps=14)

plot(plot_infections(𝕏[98], P.𝒩), plot_infections(Xtrue, P.𝒩))
























# now do mcmc, write function that also makes multipleguided_animation_unknownpar

##--------------- Next we take a really large example
# Random.seed!(58) 
# n_times = 500
# Xtrue = sample_trajectory(Ptrue::SIRforward, n_times, x0)
# pforward = plotpath(Xtrue;name="true")

# obs_times = 1:25:n_times
# 𝒪 = create_data_regular(Xtrue, obs_times, n_times, n_particles, O)



# # set guided process
# pobs = plotpath(𝒪;name="What we observe")
# ℐ =  initialise_infected_neighbours(𝒪)
# P = SIRguided(Ptrue, ℐ, 𝒪, O) 
# Π = [SA_F64[0.96, 0.04, 0.0] for _ in 1:n_particles]
# B = backward(P);


# # sample guided process 
# Z = innovations(n_times, n_particles)
# Xguided, ll  = forward(P, Π, B, Z, prior);
# ll_rounded = round(ll; digits=1)
# pguided = plotpath(Xguided; name="Reconstructed, ll is $ll_rounded")

# lo = @layout [a;b;c]
# pp = plot(pforward, pguided, pobs, layout=lo)
# png(pp,  presfigdir*"/large_example_forward_guided.png")







function mcmc_with_animation(P::SIRguided, Π, Z, prior, blocks;  δ = 0.1, γ = 0.9,
    acc = 0, ITER = 100, adaptmax=1000,
    par_estimation = false,
    propσ = 0.1,
    adaptskip = 50,
    printskip=5)
    
    # initialise
    B = backward(P)
    Zᵒ = deepcopy(Z)
    X, ll  = forward(P, Π, B, Z, prior)  # ll should also contain log prior !!!
    ℐ = P.ℐ
    
    Xᵒ = deepcopy(X)
    Bᵒ = deepcopy(B)
    XX = [copy(X)]
    ZZ = [copy(Z)]
    lls = [ll]
    θs = [param(P)]

   # mcmc
    anim = @animate  for i in 1:ITER
        ll, Z, Zᵒ, X, Xᵒ, acc_ = mcmc_iteration!(P, Π, Z, Zᵒ, X, Xᵒ, B, ll, prior, blocks, δ, i, printskip)
        acc += acc_
        acc_ = 0
        #if i ÷ 10 == 0
             push!(XX, deepcopy(X))
             push!(ZZ, deepcopy(Z))
        #end
        plotpath(X; name="")

        # adapt nr of infections in guided proposal
        if i < adaptmax #&& i÷adaptskip==0
            ℐ = γ * ℐ  + (1.0-γ) * count_infections(X, 𝒩)
        end
        if i==adaptmax
            @reset P.ℐ = ℐ
            B = backward(P)
            ll = forward!(X, P, Π, B, Z, prior)
        end

        if par_estimation
            i ÷ 10 == 0 && println(P)
            #ll, P, X, Xᵒ, B = update_parameters!(P, Xᵒ, X, Π, B, Z, ll, propσ, prior)
            ll, P, X, Xᵒ, B, Bᵒ = update_parameters!(P, X, Xᵒ, B, Bᵒ, Π, Z, ll, propσ, prior)
            push!(θs, param(P))
        end
       push!(lls, ll)
   end
   @show acc/(ITER*n_blocks)
   XX, ZZ, lls, θs, P, anim 
end






blocksize = 10
n_blocks= n_times ÷ blocksize
blocks = make_partition(n_times, n_blocks)
#blocks = [1, 2, 3:5, 6:n_times]
Xs, Zs, lls, θs, Pout, anim  = mcmc_with_animation(P, Π, Z, prior, blocks; δ=0.01,
                ITER=1500,
                γ=0.8,
                adaptmax= 00,
                adaptskip=10,
                par_estimation=false, 
                printskip=10) ;  

plot(lls)

lo = @layout [a;b;c]
plot(pforward,plotpath(Xs[end];name="Final iteration"),pobs, layout=lo)
plot_infections(Xs[end], 𝒩)

mp4(anim,presfigdir*"/mcmc_guided.mp4", fps=20)

# what happens to innovations Z?
animZ = @animate for i in eachindex(Zs)
    heatmap(hcat([a for a in Zs[i]]...))
end

mp4(animZ,presfigdir*"/mcmc_guided_Z.mp4", fps=20)

# check Cov in Z
aa=[vcat(Zs[i]...) for i in eachindex(Zs)]
aaa = hcat(aa...)'
heatmap(cov(aaa))


# with parameter par_estimation

Pinit = @set Ptrue.μ = 1.0
Pinit = @set Pinit.ν = 1.5
P = SIRguided(Pinit, ℐ, 𝒪, O) 

Xs, Zs, lls, θs, Pout, anim  = mcmc_with_animation(P, Π, Z, prior, blocks; δ=0.0001,
                ITER=2000,
                adaptmax= 100,
                adaptskip = 0,
                propσ=0.02,
                par_estimation=true, 
                printskip=25) ;  
plot(pforward,plotpath(Xs[end]),pobs, layout=lo)


pi1 = plot_infections(Xs[end], 𝒩;name="final iteration")
pi2 = plot_infections(Xtrue, 𝒩; name="true")

my_palette = ["#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"]
pi3 = heatmap(hcat([a for a in Pout.ℐ]...), color=my_palette,title="in guided process")


pi123 = plot(pi2, pi1, pi3; layout=lo, size=(500,600))
png(pi123,  presfigdir*"/reconstruction_infected.png")


bi = 1
lo = @layout [a;b;c;d]
λs = getindex.(θs,:λ);
μs = getindex.(θs,:μ);
νs = getindex.(θs,:ν);
plot(plot(lls), plot(λs[bi:end], title="λ", label=""),
                plot(μs[bi:end], title="μ", label=""),
                plot(νs[bi:end], title="ν", label=""),
                layout=lo, size=(400, 600))



# Try SMC
NUMPARTICLES = 100 
NR_SMC_STEPS = 10
NR_MOVE_STEPS =0
printskip = 1000
δstep = 0.005


blocksize = 30
n_blocks= n_times ÷ blocksize
blocks = make_partition(n_times, n_blocks)


out = smc(NR_SMC_STEPS, NUMPARTICLES, NR_MOVE_STEPS, P, Π, prior, blocks, δstep, printskip)

lo = @layout [a;b;c]
an = @animate for j in 1:NUMPARTICLES
    Xj = out.particles[j].X
    plot(pforward,plotpath(Xj),pobs, layout=lo)
end

mp4(an,presfigdir*"/smc_guided.mp4", fps=10)

plot(vcat(out.logweights...))

# try the following: have 1000 particles, don't do pcn updates
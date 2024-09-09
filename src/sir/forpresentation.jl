## make simple plots to explain the dynamics

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
    p = plotpath(Xpart)
    png(p, presfigdir*"/smallexample$i.png")
end

anim = @animate for i in 1:n_times
    Xpart = latent
    Xpart[1:i] = X[1:i]
    p = plotpath(Xpart)
    p
end
mp4(anim,presfigdir*"/smallexample_animation.mp4", fps=1)




# now a bigger example 
n_particles = 50
n_times = 100

Random.seed!(30)



# set neighbourhood structure
𝒩 = set_neighbours(n_particles, size_neighbourhood)

# set true pars
ξ, λ, μ, ν, τ =  1.0, 4.0, 2.5, 0.5, 0.1
Ptrue = SIRforward(ξ, λ, μ, ν, τ, 𝒩)

# set initial state
x0 = vcat(_I_, fill(_S_,n_particles-2),_I_) 

anim2 = @animate for i in 1:5
    Xtrue = sample_trajectory(Ptrue::SIRforward, n_times, x0)
    pforward = plotpath(Xtrue; name="forward")
    #png(pforward, presfigdir*"/large_example$i.png")
end 
mp4(anim2,presfigdir*"/largeexample_multiplerealisations_animation.mp4", fps=0.5)


# also when higher rate of getting again susceptbile
ξ, λ, μ, ν, τ =  1.0, 4.0, 2.5, 1.5, 0.1
Ptrue2 = SIRforward(ξ, λ, μ, ν, τ, 𝒩)

# set initial state
x0 = vcat(_I_, fill(_S_,n_particles-2),_I_) 

anim3 = @animate for i in 1:5
    Xtrue = sample_trajectory(Ptrue2::SIRforward, n_times, x0)
    pforward = plotpath(Xtrue; name="forward")
    #png(pforward, presfigdir*"/large_example$i.png")
end 
mp4(anim3,presfigdir*"/largeexample_multiplerealisations2_animation.mp4", fps=0.5)


## some plots of observations

# set observation scheme
δobs = 0.001
O = SA[1.0-δobs δobs/2.0 δobs/2.0; δobs/2.0 1.0- δobs δobs/2.0; δobs/2.0 δobs/2.0 1-δobs]

samplesize = (n_times * n_particles)÷20

Random.seed!(6) # nice!!
Random.seed!(666)
Xtrue = sample_trajectory(Ptrue::SIRforward, n_times, x0)
𝒪 = create_data(Xtrue, samplesize, n_times, n_particles, O)
Xobs = [O.x for O in 𝒪]


# visualise
lo = @layout [a;b]
pforward = plotpath(Xtrue;name="True, unobserved")
pobs = plotpath(Xobs;name="What we observe")
pp = plot(pforward, pobs, layout=lo)
png(pp,  presfigdir*"/large_example_forw_and_observe.png")

###############################################################
# set guided process
frac_infected_observed = sum(Xobs_flat .== _I_)/(length(Xobs_flat) - sum(Xobs_flat .== _L_))
ℐ = [fill(frac_infected_observed, n_particles) for _ ∈ 1:n_times] # of course these obs schemes use some bias but fine if only first step
#P = SIRguided(0.999, Ptrue.λ ,  Ptrue.μ, Ptrue.ν, Ptrue.τ, Ptrue.𝒩, ℐ, 𝒪, O)
P = SIRguided(1.0, Ptrue.λ ,  Ptrue.μ, Ptrue.ν, Ptrue.τ, Ptrue.𝒩, ℐ, 𝒪, O)

# sample guided process 
Random.seed!(58) 


Π = [SA_F64[0.7, 0.3, 0.0] for _ in 1:n_particles]
B = backward(P);
Z = innovations(n_times, n_particles)
Xguided, ll  = forward(P, Π, B, Z);
pguided = plotpath(Xguided; name="Reconstructed")

lo = @layout [a;b;c]
pp = plot(pforward, pguided, pobs, layout=lo)
png(pp,  presfigdir*"/large_example_forward_guided.png")

# sample multiple guided processes

Zbest = innovations(n_times, n_particles)
Xguidedbest, llbest = forward(P, Π, B, Zbest)


anim4 = @animate for i in 1:30
global llbest    
    Z = innovations(n_times, n_particles)
    Xguided, ll  = forward(P, Π, B, Z);
    lll =  round(ll; digits=1)
    pguided = plotpath(Xguided; name="Reconstructed. $lll")
    if ll > llbest
        Zbest .= Z
        llbest = ll
    end
    @show ll
    pp = plot(pforward, pguided, layout=lo)
    pp
end
@show llbest

mp4(anim4,presfigdir*"/multipleguided_animation.mp4", fps=1)

Xguidedbest, llbest = forward(P, Π, B, Zbest)
plot(pforward, plotpath(Xguidedbest), layout=lo)

# with the best path, update the guess for infected neighbours
ℐ = count_infections(Xguidedbest, 𝒩)
P2 = SIRguided(1.0, Ptrue.λ ,  Ptrue.μ, Ptrue.ν, Ptrue.τ, Ptrue.𝒩, ℐ, 𝒪, O)

B2 = backward(P2)
Xguidedbest2, llbest2 = forward(P2, Π, B2, Zbest)
plot(pforward, plotpath(Xguidedbest2), layout=lo)

plot(pforward, plotpath(Xguidedbest), plotpath(Xguidedbest2))

plot(pforward, heatmap(obs2matrix(ℐ)', title="nr infected neighbours"))









# now do the same but with parameters sampled from their prior


# sample multiple guided processes
prior = (μ=Exponential(5.0), λ = Exponential(5.0), ν=Exponential(5.0))
anim5 = @animate for i in 1:30
    λ, μ, ν = rand(prior.λ), rand(prior.μ), rand(prior.ν)
    P = SIRguided(Ptrue.ξ, λ ,  μ, ν, Ptrue.τ, Ptrue.𝒩, ℐ, 𝒪, O)
    B = backward(P)
    Z = innovations(n_times, n_particles)
    Xguided, ll  = forward(P, Π, B, Z);
    lll =  round(ll; digits=1)
    pguided = plotpath(Xguided; name="Reconstructed. $lll")

    pp = plot(pforward, pguided, layout=lo)
    pp
end
mp4(anim5,presfigdir*"/multipleguided_animation_unknownpar.mp4", fps=1)

## make simple plots to explain the dynamics
ξ = 0.999 # prob 1-ξ of new infection (not from neighbours)
τ = 0.1 # time step
δobs = 0.00

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




# now a bigger example 
n_particles = 50
n_times = 100
size_neighbourhood = 2
latent = fill(fill(_L_, n_particles),n_times)

Random.seed!(30)



# set neighbourhood structure
𝒩 = set_neighbours(n_particles, size_neighbourhood)

# set true pars
λ, μ, ν =  4.0, 2.5, 0.5
Ptrue = SIRforward(ξ, λ, μ, ν, τ, 𝒩)

# set initial state
x0 = vcat(_I_, fill(_S_,n_particles-2),_I_) 


X = sample_trajectory(Ptrue::SIRforward, n_times, x0)
anim2a = @animate for i in 1:n_times
    Xpart = latent
    Xpart[1:i] = X[1:i]
    p = plotpath(Xpart)
    p
end
mp4(anim2a,presfigdir*"/large_example.mp4", fps=10)


# other realisations
Random.seed!(2)
for i in 1:8
    Xtrue = sample_trajectory(Ptrue::SIRforward, n_times, x0)
    pforward = plotpath(Xtrue; name="")
    png(pforward, presfigdir*"/large_example$i.png")
end 



# also when higher rate of getting again susceptbile

Ptrue2 = @set Ptrue.ν = 1.5

Random.seed!(2)
for i in 1:8
    Xtrue = sample_trajectory(Ptrue2::SIRforward, n_times, x0)
    pforward = plotpath(Xtrue; name="")
    png(pforward, presfigdir*"/large_example_faster_susceptible$i.png")
end 


## some plots of observations

# set observation scheme
O = SA[1.0-δobs δobs/2.0 δobs/2.0; δobs/2.0 1.0- δobs δobs/2.0; δobs/2.0 δobs/2.0 1-δobs]



Random.seed!(666) # nice!!

# In the followoing, Omessages determines the message sent from observations
# δ = 0.01 # in the guided process assume some noise 
# Omessages = SA[1.0-δ δ/2.0 δ/2.0; δ/2.0 1.0- δ δ/2.0; δ/2.0 δ/2.0 1-δ]

Xtrue = sample_trajectory(Ptrue::SIRforward, n_times, x0)

samplesize = (n_times * n_particles)÷10   #20
𝒪 = create_data(Xtrue, samplesize, n_times, n_particles, O)

obs_times = 10:10:100
𝒪2 = create_data_regular(Xtrue, obs_times, n_times, n_particles, O)

obs_particles = 5:7:n_particles 
𝒪3 = create_data_see_some_particles(Xtrue, obs_particles, n_times, n_particles, O)

#


# visualise
lo = @layout [a;b]
pforward = plotpath(Xtrue;name="True, unobserved")
pobs = plotpath(𝒪;name="What we observe")
pp = plot(pforward, pobs, layout=lo)
png(pp,  presfigdir*"/large_example_forw_and_observe1.png")

pobs2 = plotpath(𝒪2;name="What we observe")
pp = plot(pforward, pobs2, layout=lo)
png(pp,  presfigdir*"/large_example_forw_and_observe2.png")



###############################################################
Random.seed!(58) 

# set guided process
𝒪 = 𝒪2
pobs = plotpath(𝒪;name="What we observe")
ℐ =  initialise_infected_neighbours(𝒪)
P = SIRguided(Ptrue, ℐ, 𝒪, O) 
Π = [SA_F64[0.96, 0.04, 0.0] for _ in 1:n_particles]
B = backward(P);


# sample guided process 
Z = innovations(n_times, n_particles)
Xguided, ll  = forward(P, Π, B, Z);
ll_rounded = round(ll; digits=1)
pguided = plotpath(Xguided; name="Reconstructed, ll is $ll_rounded")

lo = @layout [a;b;c]
pp = plot(pforward, pguided, pobs, layout=lo)
png(pp,  presfigdir*"/large_example_forward_guided.png")

# sample multiple guided processes

Zbest = innovations(n_times, n_particles)
Xguidedbest, llbest = forward(P, Π, B, Zbest)

lo = @layout [a;b]
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
    pp = plot(pforward, pguided, layout=lo)#, size=(400, 790))
    pp
end
@show llbest

mp4(anim4,presfigdir*"/multipleguided.mp4", fps=2)


Xguidedbest, llbest = forward(P, Π, B, Zbest)
lo = @layout [a;b;c]
plot(pforward, plotpath(Xguidedbest), pobs, layout=lo)

# now do mcmc, write function that also makes multipleguided_animation_unknownpar

function mcmc_with_animation(P::SIRguided, Π, Z, blocks;  δ = 0.1, γ = 0.9,
    acc = 0, ITER = 100, adaptmax=1000,
    par_estimation = false,
    propσ = 0.1,
    prior = (μ=Exponential(5.0), λ = Exponential(5.0), ν=Exponential(5.0)),
    adaptskip = 50,
    printskip=5)
    

            #    @unpack 𝒪 = P
                #    n_times, n_particles = length(𝒪), length(𝒪[1].x)
            #n_blocks = length(blocks)

   B = backward(P)
   Zᵒ = deepcopy(Z)
   X, ll  = forward(P, Π, B, Z)
   Xᵒ = deepcopy(X)

   XX = [copy(X)]
   lls = [ll]
   θs = [param(P)]

    anim = @animate  for i in 1:ITER
       for block in blocks
           update!(Zᵒ, Z, δ, block)
           llᵒ = forward!(Xᵒ, P, Π, B, Zᵒ)
           if log(rand()) < llᵒ - ll
               mod(i,printskip)==0 && println(i,"  ",ll,"  ", llᵒ,"  ", llᵒ-ll, "  accepted")
               ll = llᵒ
               Z[block] .= Zᵒ[block]
               X .= Xᵒ
               acc += 1
           else
               mod(i,printskip)==0 && println(i, "   ", ll,"  ", llᵒ,"  ", llᵒ-ll, "  rejected")
           end
        end
        i ÷ 10 == 0 && push!(XX, deepcopy(X))
        #i ÷ 10 == 0 && plotpath(X; name="")
        plotpath(X; name="")
       

       if i < adaptmax && i÷adaptskip==0
           ℐnew = γ * P.ℐ  + (1.0-γ) * count_infections(X, 𝒩)
           @reset P.ℐ = ℐnew
           B = backward(P)
           ll = forward!(X, P, Π, B, Z)
       end

       if par_estimation
           # update μ
           μᵒ = P.μ * exp(propσ*randn())
           logprior_proposalratios = (log(μᵒ) - log(P.μ)) + logpdf(prior.μ,μᵒ) - logpdf(prior.μ,P.μ)
           Pᵒ = @set P.μ = μᵒ
           ll, P = updatepar!(Xᵒ, X, Pᵒ, P, Π,  B, Z,  ll, logprior_proposalratios)

           # update λ
           λᵒ = P.λ * exp(propσ*randn())
           logprior_proposalratios = (log(λᵒ) - log(P.λ)) + logpdf(prior.λ,λᵒ) - logpdf(prior.λ,P.λ)
           Pᵒ = @set P.λ = λᵒ
           ll, P = updatepar!(Xᵒ, X, Pᵒ, P, Π,  B, Z,  ll, logprior_proposalratios)


           # update ν
           νᵒ = P.ν * exp(propσ*randn())
           logprior_proposalratios = (log(νᵒ) - log(P.ν)) + logpdf(prior.ν,νᵒ) - logpdf(prior.ν,P.ν)
           Pᵒ = @set P.ν = νᵒ
           ll, P = updatepar!(Xᵒ, X, Pᵒ, P, Π,  B, Z, ll, logprior_proposalratios)
           push!(θs, param(P))
       end
       push!(lls, ll)
   end
   @show acc/(ITER*n_blocks)
   XX, lls, θs, P, anim 
end






blocksize = 10
n_blocks= n_times ÷ blocksize
blocks = make_partition(n_times, n_blocks)
#blocks = [1, 2, 3:5, 6:n_times]
Xs, lls, θs, Pout, anim  = mcmc_with_animation(P, Π, Zbest, blocks; δ=0.005,
                ITER=200,
                γ=0.8,
                adaptmax= 200,
                adaptskip=10,
                par_estimation=false, 
                printskip=1) ;  



plot(pforward,plotpath(Xs[end];name="Final iteration"),pobs, layout=lo)
plot_infections(Xs[end], 𝒩)

mp4(anim,presfigdir*"/mcmc_guided.mp4", fps=20)

# with parameter par_estimation

Pinit = @set Ptrue.λ = 2.0
P = SIRguided(Pinit, ℐ, 𝒪, O) 



Xs, lls, θs, Pout, anim  = mcmc_with_animation(P, Π, Zbest, blocks; δ=0.005,
                ITER=200,
                adaptmax= 1000,
                adaptskip = 20,
                propσ=0.02,
                par_estimation=true, 
                printskip=5) ;  
plot(pforward,plotpath(Xs[end]),pobs, layout=lo)

lo = @layout [a;b;c;d]
λs = getindex.(θs,:λ);
μs = getindex.(θs,:μ);
νs = getindex.(θs,:ν);
plot(plot(lls), plot(λs, title="λ", label=""),
                plot(μs, title="μ", label=""),
                plot(νs, title="ν", label=""),
                layout=lo, size=(400, 600))










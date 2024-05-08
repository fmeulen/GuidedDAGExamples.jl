using StaticArrays
abstract type MarkovProcess end

struct SIRforward <: MarkovProcess
        ξ::Float64
        λ::Float64
        μ::Float64
        ν::Float64
        τ::Float64
        𝒩::Array{Array{Int64,1},1}
end

struct SIRguided{T} <: MarkovProcess
    ξ::Float64
    λ::Float64
    μ::Float64
    ν::Float64
    τ::Float64
    𝒩::Array{Array{Int64,1},1}
    ℐ::Array{Array{T,1},1}  # vector with each element a vector at that time of nr of infected neighbours
end

@enum State::UInt8 _S_=1 _I_=2 _R_=3 _L_=0 
const 𝒳 = @SVector [_S_,_I_,_R_]

# state space of one person
const BF = [@SVector([3,1,2]), @SVector([1,2,3]), @SVector([2,3,1])]


######## function defns

params(P::MarkovProcess) = (λ=P.λ, μ=P.μ, ν=P.ν)

"""
    ind(x)

x ∈ 𝒳 is mapped into integer
"""
ind(x) = Int(x)


function observationmessage(x::State)
    if x==_S_ 
        return(SA_F64[1, 0, 0])
    elseif x==_I_
        return(SA_F64[0, 1, 0])
    elseif x==_R_
        return(SA_F64[0, 0, 1])
    else 
        return(SA_F64[1, 1, 1])  
    end
end

"""
    set_neighbours(n)

Returns an array of length n, where the i-th element contains the indices for the neighbours of individual i
"""
function set_neighbours(n, k)
    if k==2
        if n<5
        @error "provide larger value for n"
        end
        𝒩 = Vector{Vector{Int64}}(undef,n)
        𝒩[1] = [2,3]; 𝒩[2] = [1,3,4]; 𝒩[n-1]  = [n, n-2, n-3]; 𝒩[n] = [n-1,n-2]
        for i in 3:n-2
            𝒩[i] = [i-2,i-1,i+1,i+2]
        end
    end
    if k==1
        𝒩 = Vector{Vector{Int64}}(undef,n)
        𝒩[1] = [2]; 𝒩[n] =[n-1]
        for i in 2:n-1
            𝒩[i] = [i-1, i+1]
        end
    end
    𝒩
end

#### forward sampling #####
rand𝒳(z, p)= z < p[1] ? _S_ : ( z > 1-p[3] ? _R_ : _I_ )

"""
    sample𝒳(x,z,p)

Sample from 𝒳 according to probability vector p

p =  [1/4, 1/2, 1/4]
z = rand()
sample𝒳(_S_, z, p)
"""
function sample𝒳(x, u::Float64,p) # provide current state
    b, f = BF[Int(x)][1], BF[Int(x)][3]
    if u < p[b] return(𝒳[b])
    elseif u > 1-p[f] return(𝒳[f])
    else return x
    end
end

pS(x) = @SVector [exp(-x), 1.0-exp(-x), 0.0]  # λ*τ
pI(x) = @SVector [0.0, exp(-x), 1.0-exp(-x)]      # μ*τ
pR(x) = @SVector [1.0-exp(-x), 0.0, exp(-x)]     # ν*τ



function κ(P::SIRforward, i, x) # used to be called pi
    if x[i] == _S_
        pS(P.λ * nr_infected_neighb(x,P.𝒩,i) * P.τ)
    elseif x[i] == _I_
        pI(P.μ * P.τ)
    elseif x[i] == _R_
        pR(P.ν * P.τ)
    end
end

"""
    sample_particle(P::SIRforward,i,x,z)

One forward simulation step for invdividual `i`, if the present configuration of all individuals is `x`
i: index of individual to be forward simulated
x: current state of all individuals
z: innovation for this step
"""
function sample_particle(P::SIRforward, i, x, z) 
    p = κ(P, i, x)
    sample𝒳(x[i], z, p)
end

sample(P::SIRforward, x, z) = [sample_particle(P, i, x, z[i]) for i in eachindex(x)]

"""
    sample_trajectory(P::SIRforward, n::Int, x0)

    sample SIR-process over n time instances, where the first time instance is x0
"""
function sample_trajectory(P::SIRforward, n_times::Int64, x0)
    X = [x0]
    n_particles = length(x0)
    for j in 2:n_times
        z = rand(n_particles)
        push!(X, sample(P, X[j-1], z))
    end
    X
end







obs2matrix(X) =  [ind(X[j][i]) for j in eachindex(X), i in eachindex(X[1])]


function updatepars!(P, Pᵒ ,X ,Xᵒ ,lr ,lrᵒ, Xobs, Z, Q ,propσ, prior, accpar, it, skip_print)
    nseg = length(Xobs)-1
    λᵒ = P.λ * exp(propσ*randn())
    Pᵒ = SIRguided(λᵒ, P.μ, P.ν, P.τ, P.𝒩)
    for k in 1:nseg
        Xᵒ[k], lrᵒ[k] = sample_segment!(Pᵒ, Xᵒ[k], Xobs[k],Xobs[k+1],Z[k],Q[k],J)
    end
    if log(rand()) < sum(lrᵒ .- lr)  + (log(Pᵒ.λ) - log(P.λ)) + logpdf(prior[1],λᵒ) - logpdf(prior[1],P.λ)
        if mod(it, skip_print)==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    λ ", λᵒ)  end
        X, Xᵒ, lr, lrᵒ, P, Pᵒ = Xᵒ, X, lrᵒ, lr, Pᵒ, P
        accpar[1] += 1
    else
        if mod(it, skip_print )==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    λ ", P.λ)  end
    end

    μᵒ = P.μ  * exp(propσ*randn())
    Pᵒ = SIRguided(P.λ, μᵒ, P.ν, P.τ, P.𝒩)
    for k in 1:nseg
        Xᵒ[k], lrᵒ[k] = sample_segment!(Pᵒ, Xᵒ[k], Xobs[k],Xobs[k+1],Z[k],Q[k],J)
    end
    if log(rand()) < sum(lrᵒ .- lr)  + (log(Pᵒ.μ) - log(P.μ))+ logpdf(prior[2],μᵒ) - logpdf(prior[2],P.μ)
        if mod(it, skip_print)==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    μ ", μᵒ)  end
        X, Xᵒ, lr, lrᵒ, P, Pᵒ = Xᵒ, X, lrᵒ, lr, Pᵒ, P
        accpar[2] += 1
    else
        if mod(it, skip_print )==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    μ ", P.μ)  end
    end

    νᵒ = P.ν  * exp(propσ*randn())
    Pᵒ = SIRguided(P.λ, P.μ, νᵒ, P.τ, P.𝒩)
    for k in 1:nseg
        Xᵒ[k], lrᵒ[k] = sample_segment!(Pᵒ, Xᵒ[k], Xobs[k],Xobs[k+1],Z[k],Q[k],J)
    end
    if log(rand()) < sum(lrᵒ .- lr)  + (log(Pᵒ.ν) - log(P.ν)) + logpdf(prior[3],νᵒ) - logpdf(prior[3],P.ν)
        if mod(it, skip_print)==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    ν ", νᵒ)  end
        X, Xᵒ, lr, lrᵒ, P, Pᵒ = Xᵒ, X, lrᵒ, lr, Pᵒ, P
        accpar[3] += 1
    else
        if mod(it, skip_print )==0 println("iteration ", it,  "   diff loglr ", sum(lrᵒ.-lr), "    ν ", P.ν)  end
    end


    return X, Xᵒ, lr, lrᵒ, P, Pᵒ

end

# convention:
# k indexes segments k ∈ 1...n_times-1
# j indexes time on a segment j ∈ 1...J
# i indexes person i ∈ 1...n


function sir_inference(Xobs, P, J; ρ = 0.99, propσ = 0.1,
    ITER = 100, skip_print = 10,
    prior =(Exponential(1.0), Exponential(1.0), Exponential(1.0)),
    γ = 0.5, adaptmax = div(1*ITER,3))

    θ = params(P)
    θθ = [θ]
    n_times = length(Xobs)
    n_particles = length(Xobs[1])

    X = [ [[_S_ for  i in 1:n] for j in 1:J] for _ in 2:n_times]

    # ninfected = 2mean(obs2matrix(Xobs).==ind(_I_))
    # N = [ [[ninfected for  i in 1:n] for j in 1:J] for k in 2:n_times]


    N = [ [[float.(nr_infected_neighb(Xobs[k],P.𝒩,i)) for  i in 1:n] for j in 1:J] for k in 2:n_times]
    Q = [ [[SMatrix{3,3}(1.0I) for  i in 1:n] for j in 1:J] for _ in 2:n_times]
    h̃!(Q,θ,N,P.τ)

    lr = zeros(n_times-1)
    Z = [[rand(n) for _ in 1:J-1] for _ in 1:n_times]
    for k in 2:n_times
        _, lr[k-1] = sample_segment!(P, X[k-1], Xobs[k-1], Xobs[k], Z[k-1], Q[k-1], J)
    end


    Xinit = deepcopy(X) # save
    Xmid = 0            # save middle iteration
    Xᵒ = deepcopy(X)
    Zᵒ = deepcopy(Z)
    Pᵒ = deepcopy(P)
    lrᵒ = deepcopy(lr)

    acc = 0
    accpar = fill(0,3)
    difflr = Float64[]

    for it in 2:ITER
        X, Xᵒ, lr, lrᵒ, P, Pᵒ = updatepars!(P, Pᵒ, X, Xᵒ, lr, lrᵒ, Xobs,Z,Q,propσ,prior,accpar,it,skip_print)
        θ = params(P)
        push!(θθ, θ)

        # recompute likelihood under the auxiliary process with updated θ
        h̃!(Q, θ, N, P.τ)
        for k in 1:n_times-1
            X[k], lr[k] = sample_segment!(P, X[k], Xobs[k], Xobs[k+1], Z[k], Q[k],J)
        end

        # update innovations Z
        for k in 1:n_times-1
            Znew = [randn(n) for _ in 1:J-1]
            Zᵒ[k] = [mod.(Z[k][j] + (1.0-ρ)*Znew[j], 1) for j in eachindex(Z[1])]
            Xᵒ[k], lrᵒ[k] = sample_segment!(P, Xᵒ[k], Xobs[k], Xobs[k+1], Zᵒ[k], Q[k], J)
            Δlr = lrᵒ[k] - lr[k]
            push!(difflr, Δlr)
            if log(rand()) < Δlr
                acc += 1
                X[k] = Xᵒ[k]
                Z[k] = Zᵒ[k]
                lr[k] = lrᵒ[k]
            end
            if (mod(it, skip_print)==0) & (k==1)
                println("iteration ", it   ,"   diff loglr ", Δlr)
                println("----------")
             end
        end

        #    update N
        # if adaptmin < it < adaptmax
        #     for k in 1:n_times-1
        #         N_ =  [[adaptfrac*ninfected + (1-adaptfrac)*nr_infected_neighb(X[k][j], P.𝒩, i) for  i in 1:n] for j in 1:J]
        #         N[k] = N_/it + (it/(it+1))*  N[k]
        #     end
        #     h̃!(Q,θ,N,P.τ)
        # end
        if it < adaptmax

                N_ = [  [[float.(nr_infected_neighb(X[k][j], P.𝒩, i)) for  i in 1:n] for j in 1:J] for k in 1:n_times-1]
                #N = [ [[float.(nr_infected_neighb(Xobs[k],P.𝒩,i)) for  i in 1:n] for j in 1:J] for k in 2:n_times]
                N .= γ * N + (1-γ) * N_

            h̃!(Q,θ,N,P.τ)
        end

        if it==div(ITER,2)
            Xmid = deepcopy(X)
        end
#        println(Xmid)
    end
    θθ, X, N, Q, Xinit, Xmid, sum(acc)/(n_times*ITER), sum.(accpar)/ITER, difflr
end

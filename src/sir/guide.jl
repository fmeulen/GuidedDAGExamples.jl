
##########################
# B is backward filter output
function guide(x, P::SIRguided, h, z)
    @assert length(x)==length(h)==length(z)
    y = similar(x)
    logweight = 0.0 
    for i ∈ eachindex(x)
        if x[i]==_S_
            ni = nr_infected_neighb(x, P.𝒩, i)
            p = pS(P.λ * ni * P.τ) .* h[i]
        elseif x[i]==_I_
            p = pI(P.μ * P.τ) .* h[i]
        elseif x[i]==_R_
            p = pR(P.ν * P.τ) .* h[i]
        end
        y[i] = rand𝒳(z[i], p) #sample𝒳(x[i], z[i], p/sum(p))
        logweight += log(sum(p))
    end
    y, logweight
end

rand𝒳(z, p)= z < p[1] ? _S_ : ( z > 1-p[3] ? _R_ : _I_ )
    

function forward(P::SIRguided, Π, B, Z)
    # sample initial state
    X = Vector{State}(undef, n_particles)
    z = Z[1]
    for i in 1:n_particles
        p = Π[i] .* B[1][i]
        X[i] = rand𝒳(z[i], p/sum(p))
    end
    ll = log(sum(p))

    Xs = [X]
    for t in 2:n_steps
        Xnext, lw = guide(X, P, B[t], Z[t])
        X = Xnext
        ll += lw
        push!(Xs, copy(X))
    end
    Xs, ll
end


"""
one step transition prob for P
"""
function logP(P::MarkovProcess, x,y)
    out = 0.0
    for i in eachindex(x)
        iy = ind(y[i])
        if x[i]==_S_
            ni = nr_infected_neighb(x,P.𝒩,i)
            out += log(pS(P.λ * ni * P.τ)[iy])
        elseif x[i]==_I_
            out+= log(pI(P.μ * P.τ)[iy])
        elseif x[i]==_R_
            out+= log(pR(P.ν * P.τ)[iy])
        end
    end
    out
end



function sample_segment!(P::SIRguided, Xᵒ, xstart, xend, Zseg, Qseg, J)
    Xᵒ[1] = xstart
    logLR = 0.0
    for j in 2:J
        xᵒ, logwj = Pstep(P, Xᵒ[j-1], Zseg[j-1], Qseg[j], xend)
        Xᵒ[j] = xᵒ
        logLR += logwj
    end
    logLR += logP(P, Xᵒ[J-1], Xᵒ[J])
    Xᵒ, logLR
end

################ simulating guided proposal ####################################

"""
    Q̃(θ,ninfected,τ)

Make Q̃ matrix with parameter θ which is assumed a named tuple with elements λ, μ and ν;
ninfected is the number of infected neighbours, τ is the time-discretisation step
"""

Q̃(θ,ninfected,τ) = hcat(pS(θ.λ * τ * ninfected), pI(θ.μ*τ), pR(θ.ν*τ))'

"""
    Qᵒstep(P::SIRguided,i,x,zi, q̃)

Returns state at time j for individual i as well as its log-weight
"""
function Qstep(P::SIRguided,i,x,zi, q̃)
    out =_R_; p = 0.0
    if x[i]==_S_
        ni = nr_infected_neighb(x, P.𝒩, i)
        p = pS(P.λ * ni * P.τ) .* q̃
    elseif x[i]==_I_
        p = pI(P.μ * P.τ) .* q̃
    elseif x[i]==_R_
        p = pR(P.ν * P.τ) .* q̃
    end
    out = sample𝒳(x[i], zi, p/sum(p))
    out, log(sum(p))
end

"""
    Pstep(x,θ,𝒩,xend,τ,θ̃,j,J)

step from time j-1 to j, all individuals
"""
function Pstep(P::SIRguided,x,z,Qseg,xend) # Qseg contains all matrices for that segment
    xᵒ = State[]
    logw = 0.0
    for i in eachindex(x)
        Q̃mat = Qseg[i]   #      #Q̃Jj = Q̃(θ̃,τ)^(J-j)      Q̃Jj = prod([Q̃((λ=Pᵒ.λ*Nj[k][i], μ=Pᵒ.μ, ν=Pᵒ.ν),Pᵒ.τ) for k in (j+1):J])
        i_endstate = ind(xend[i])
        xnext, w = Qstep(P,i,x,z[i],Q̃mat[:,i_endstate])
        push!(xᵒ, xnext)
        logw += w
        logw -= log(Q̃mat[ind(xnext), i_endstate])
    end
    xᵒ, logw
end

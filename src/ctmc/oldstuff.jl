
########### ignore everyting below here ###################

if false

# new try
#ζ(x, R::GuidedReactionNetwork) = [LinearAlgebra.dot(R.xT-x,R.ξ[ℓ]) for ℓ ∈ eachindex(R.ξ)]
ζ(s, x, R::GuidedReactionNetwork) = [LinearAlgebra.dot(R.xT - x - (R.T-s)*R.μ, R.ξ[ℓ]) /λ(R.xT, R)[ℓ] for ℓ ∈ eachindex(R.ξ)]
λ(s, x, R::GuidedReactionNetwork) = λ(x, R) .* ζ(s, x, R) /(R.T - s)

function simstep!(T, P, x, t, R::GuidedReactionNetwork)
    z = ζ(t, x, R)
    imin = findall(x-> x < 0, z)
    ipos = findall(x-> x >= 0, z)
    a = λ(x,R)
    #aᵒ = λ(t,x,R)
    Δt = abs.((P - T) ./ a)
    for ℓ ∈ ipos
        Δt[ℓ] = (R.T-t) * (1-exp(-Δt[ℓ]/z[ℓ]))
    end
    (Δ, μ) = findmin(Δt)
    for ℓ ∈ imin
        T[ℓ] += Δ * a[ℓ]  # update internal times ## ADJUST
    end
    for ℓ ∈ ipos
        T[ℓ] += a[ℓ] * z[ℓ] *log((R.T-t)/(R.T-t-Δ))
    end
    t += Δ  # update clock time
    x += R.ξ[μ]   # update state according to reaction taking place
    P[μ] += randexp()
    T, P, x, t
end

end
function τleap(t0,x0, R::GuidedReactionNetwork; h=0.1)
    t = t0
    x = x0
    eventtimes = [t]
    eventvals = [x]
    while t < R.T
        L = λ(t,x,R)
    #    println(L)
        for ℓ ∈ 1:R.nreact
            x += R.ξ[ℓ] * Random.rand(Poisson(h*L[ℓ]))
        end
        t += h
        push!(eventtimes, t)
        push!(eventvals, x)
    end
    eventtimes, eventvals
end


# some old stuff.
if false
    # ei(x) = expint(x)[1]
    # ei1(x)= expint_E1(x)[1]
    # g(ℓ, u, α) = α>0 ? u*exp(-α/u) - ℓ*exp(-α/ℓ) + α*(ei(α/u) - ei(α/ℓ)) : u*exp(-α/u) - ℓ*exp(-α/ℓ) - α*(ei1(-α/ℓ) - ei1(-α/u))
    G(λℓ, ψℓ, Tmint, PℓminTℓ) = (Δ) -> λℓ * g(Tmint-Δ, Tmint, -ψℓ) - PℓminTℓ
    ∇G(λℓ, ψℓ, t,R::GuidedReactionNetwork) = (Δ) -> -λℓ *   exp(ψℓ/(R.T-t-Δ))
end

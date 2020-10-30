Ei(x) = -real(expint(Complex(-x)))
g(u, a) = exp(-a/u)*u + a*Ei(-a/u)

"""
    g(l, u, a)

Integrate  exp(-a/x) from l to u
"""
g(l, u, a) = g(u, a) - g(l, a)


g(1, 2, pi) # should be 0.124
# https://www.wolframalpha.com/input/?i=integrate+exp%28-pi%2Fu%29+from+1+to+2


function rootsolve(λ_, s, ψ_, c)
    U = H(λ_, s, ψ_, c)
    if U(0.0) > 0
        z = Inf
    elseif U(s-0.001) < 0
        z = Inf
    else
        # ∇U = ∇H(λ_, s, ψ_, c)
        # z = find_zero((U,∇U), 0.0, Roots.Newton())
        z= find_zero(U, 0.5*s, Order1())
    end
    z
end



H(λ, s, ψ, c) = (h) -> λ * g(s-h, s, -ψ) - c
∇H(λ, s, ψ, c) = (h) -> λ * exp(ψ/(s-h))

λ_ = 1.0
s = 0.5
ψ_ = -2.0
c = 1.0
U = H(λ_, s, ψ_, c)
∇U = ∇H(λ_, s, ψ_, c)
U(0.04)

hvals = 0.01:0.01:(s-0.1)
Hvals = U.(hvals)
plot(hvals, Hvals)
∇U(0.001)
U(0.001)
rootsolve(λ_, s, ψ_, c)

z = find_zero((U,∇U), 0.0, Roots.Newton())
U(z)




rootsolve(12.0, 2.0, -ψ_, 1.0)

fzero(U,0.0)

fzero(U,0.0)
#plot(H(λ_, s, ψ_, c), 0, s)



# some old stuff.
if false
    # ei(x) = expint(x)[1]
    # ei1(x)= expint_E1(x)[1]
    # g(ℓ, u, α) = α>0 ? u*exp(-α/u) - ℓ*exp(-α/ℓ) + α*(ei(α/u) - ei(α/ℓ)) : u*exp(-α/u) - ℓ*exp(-α/ℓ) - α*(ei1(-α/ℓ) - ei1(-α/u))
    G(λℓ, ψℓ, Tmint, PℓminTℓ) = (Δ) -> λℓ * g(Tmint-Δ, Tmint, -ψℓ) - PℓminTℓ
    ∇G(λℓ, ψℓ, t,R::GuidedReactionNetwork) = (Δ) -> -λℓ *   exp(ψℓ/(R.T-t-Δ))
end

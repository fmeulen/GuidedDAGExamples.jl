Ei(x) = -real(expint(Complex(-x)))
g(u, a) = exp(-a/u)*u + a*Ei(-a/u)

"""
    g(l, u, a)

Integrate  exp(-a/x) from l to u
"""
g(l, u, a) = g(u, a) - g(l, a)


g(1, 2, pi) # should be 0.124
# https://www.wolframalpha.com/input/?i=integrate+exp%28-pi%2Fu%29+from+1+to+2
H(λ, s, ψ, c) = (h) -> λ * g(s-h, s, -ψ) - c


λ_ = 1.0
s = 0.5
ψ_ = 2.0
c = 1.0
H(λ_, s, ψ_, c)(0.04)

hvals = 0.01:0.01:(s-0.01)
H(λ_, s, ψ_, c).(hvals)

plot(H(λ_, s, ψ_, c), 0, s)

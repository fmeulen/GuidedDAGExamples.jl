
using Distributions
using Makie
using Random
g(x, α, β) = pdf(Gamma(α, 1/β), x)
lg(x, α, β) = logpdf(Gamma(α, 1/β), x)

Random.seed!(1)

T = 100
n = 1000
S = 1:n
α = (T/n)ones(n)
β = [0.1 + sin(x) for x in pi*S/n]

x = x0 = 0.0
t = 0.0
xs = [x0]
ts = [t]
for s in S
    global t, x
    t += α[s]
    x = x + rand(Gamma(α[s], 1/β[s]))
    push!(ts, t)
    push!(xs, x)

end

# save endpoint
xv = x[end]


# Gamma process bridge

x = x0 = 0.0
t = 0.0
xs2 = [x0]
ll = 0.0
βt = 0.4
for s in S[1:end-1]
    global t, x, ll
    t += α[s]
    dx = (xv - x)*rand(Beta(α[s], T - t))
    x = x + dx
    ll += lg(dx, α[s], β[s]) -  lg(dx, α[s], βt)
    push!(xs2, x)
end
push!(xs2, xv)
println("ll = $ll")

p1 = lines(ts, xs2);
lines!(p1, ts, xs, color=:red);

"""
Sample scaled exponentially shifted Beta using
rejection sampling.
"""
function randexpbeta(κ, t1, t2, s, M)
    while true
        x = s*rand(Beta(t1, t2))
        p = exp(-κ*x)/M
        @assert p <= 1
        rand() < p && return x
    end
end

function Iest(t1, t2, κ, M)
    m = 0.0
    for i in 1:M
        x = rand(Beta(t1, t2))
        m += exp(-κ*x)
    end
    return m/M
end


x = x0 = 0.0
t = 0.0
xs3 = [x0]
βt = 0.3
ll2 = 0.0
for s in S[1:end-1]
    global t, x, ll2
    dx = randexpbeta(β[s] - βt, α[s], T - t - α[s], xv - x, 1000)
    ll2 += α[s]*log(β[s]/βt) + lg(xv - x, T - t, βt) - lg(xv - (x + dx), T- t - α[s], βt)
    ll2 += log(Iest(α[s], T - t - α[s], (β[s] - βt)*(xv-x), 100))
    t += α[s]
    x = x + dx

    push!(xs3, x)
end
println("")
println("ll2 = $ll2")

push!(xs3, xv)
lines!(p1, ts, xs3, color=:blue);
p1
save("gammaprop.png", p1)
println("")

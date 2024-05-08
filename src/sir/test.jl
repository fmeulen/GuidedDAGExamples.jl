

function f(x)
    push!(x, 4.0)
    y = sin.(x)
    y
end

x = [1.0, 2.0]
f(x)
x

using ConcreteStructs

@concrete struct  Mystruct
    a
end

UU= Mystruct(Z)
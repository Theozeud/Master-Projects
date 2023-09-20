include("main.jl")

# Not forgot to declare it as constants
L = 1
Nₓ = 300
N = 10
Nplot = 100
X = range(-L/2,L/2,Nplot)

f₁(x::Real,::Real) = cos(2π*x/L)
f₂(x::Real,::Real) = exp(f₁(x,L))

f̃₁= computeApprox(f₁,L,Nₓ,N)

plotψ(f̃₁,f₁,X)

f̃₂= computeApprox(f₂,L,Nₓ,N)
plotψ(f̃₂,f₂,X)

# the rate between N and Nₓ must be enough small, i.e Nₓ>>N to have a good result
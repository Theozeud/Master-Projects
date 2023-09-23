include("discretization.jl")
include("fourier.jl")
include("plot.jl")

# Not forgot to declare it as constants
# Parameters for simulation
L = 1
Nₓ = 100
Δx = L/Nₓ
T = 10
Nₜ = 10
Δt = T/Nₜ
N = 50

# Discretisation
X = uniformDiscretization(-L/2,L/2,Nₓ)

# Somes functions for testing the code
f₁(x::Real,::Real) = cos(2π*x/L)
f₂(x::Real,t::Real) = exp(f₁(x,t))
∂f₁(x::Real,::Real) = -2π/L * sin(2π*x/L)
∂f₂(x::Real,t::Real) = ∂f₁(x,t) * exp(f₁(x,t))
∂∂f₁(x::Real,::Real) = -(2π/L)^2 * cos(2π*x/L)
∂∂f₂(x::Real,t::Real) = (∂∂f₁(x,t) + ∂f₁(x,t) )* exp(f₁(x,t))

# Test for the computation of the approximation
f̃₁= computeApprox(f₁,X,N)
@time plotψ(f₁,f̃₁,X)

f̃₂= computeApprox(f₂,X,N)
plotψ(f₂,f̃₂,X)

# Test for the computation of the approximation of the derivative
∂f̃₁= computeApprox∂(2,f₁,X,N)
plotψ(∂∂f₁,∂f̃₁,X)

∂f̃₂= computeApprox∂(2,f₂,X,N)
plotψ(∂∂f₂,∂f̃₂,X)

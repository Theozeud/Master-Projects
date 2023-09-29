include("discretization.jl")
include("fourier.jl")
include("plot.jl")

# Parameters for simulation
L = 1
Nₓ = 100
Δx = L/Nₓ
T = 10
Nₜ = 50
Δt = T/Nₜ
N = 1

# Discretisation
X = uniformDiscretization(-L/2,L/2,Nₓ)

# Somes functions for testing the code
f₁(x::Real,::Real) = cos(2π*x/L)
f₂(x::Real,t::Real) = exp(f₁(x,t))
∂f₁(x::Real,::Real) = -2π/L * sin(2π*x/L)
∂f₂(x::Real,t::Real) = ∂f₁(x,t) * exp(f₁(x,t))
∂²f₁(x::Real,::Real) = -(2π/L)^2 * cos(2π*x/L)
∂²f₂(x::Real,t::Real) = (∂²f₁(x,t) + ∂f₁(x,t)^2)* exp(f₁(x,t))

# Test for the computation of the approximation

#f̃₁= approx(f₁,X,N)
#plotψ(f₁,f̃₁,X)

#f̃₂= approx(f₂,X,N)
#plotψ(f₂,f̃₂,X)

#fftf̃₁= FFTapprox(f₁,X,1.0)
#plotψ(f₁,fftf̃₁[1],X)

#fftf̃₂= FFTapprox(f₂,X,1.0)
#plotψ(f₂,fftf̃₂[1],X)

# Test for the computation of the approximation of the derivative
#∂f̃₁= approx∂(1,f₁,X,N)
#plotψ(∂f₁,∂f̃₁,X)

#∂f̃₂= approx∂(1,f₂,X,N)
#plotψ(∂f₂,∂f̃₂,X)

#∂²f̃₁= approx∂(2,f₁,X,N)
#plotψ(∂²f₁,∂²f̃₁,X)

#∂²f̃₂= approx∂(2,f₂,X,N)
#plotψ(∂²f₂,∂²f̃₂,X)

#∂²fftf̃₁= FFTapprox∂(2, f₁,X,1.0)
#plotψ(∂²f₁,∂²fftf̃₁[1],X)

∂²fftf̃₂= FFTapprox∂(2, f₂,X,1.0)
plotψ(∂²f₂,∂²fftf̃₂[1],X)

#k = Kmatrix(0, N, X)
#plotψ(f₁,k*f₁.(X.points,0),X)


# Plots for cₖ


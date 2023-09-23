using FFTW
using SparseArrays

# Base of L²
eₖ(k::Int, x::Real, L::Real) = exp(eltype(x)(2π) * im * k *x /L )

# Function to compute the Fourier coefficent of ψ : cₖ = <eₖ,ψ>
cₖ(ψ::Base.Callable, t::Real, k::Int, X::Discretization) = 1/length(X) * eₖ.(k,-X.points,X.sup-X.inf)'*ψ.(X.points,t)

# Function to compute the Fourier coefficent of any derivative of ψ : derivcₖ = <eₖ,∂ψ>
derivcₖ(order::Int, ψ::Base.Callable, t::Real, k::Int, X::Discretization) = (eltype(X)(2π)*im*k/(X.sup-X.inf))^order * cₖ(ψ, t, k, X)

computeApprox(ψ::Base.Callable, X::Discretization, N::Int) = (x,t) -> mapreduce(k->eₖ(k,x,X.sup-X.inf)*cₖ(ψ, t, k, X),+, range(-N,N))

computeApprox∂(order::Int, ψ::Base.Callable, X::Discretization, N::Int) = (x,t) -> mapreduce(k->eₖ(k,x,X.sup-X.inf)*derivcₖ(order,ψ, t, k, X),+, range(-N,N))

# Some operators




# Functions to solve the schrodinegr equation

expmatrix(K::Matrix, ψ₀::Vector, t::Real) = exp(-im*K*t)*ψ₀

solveschrodinger(ψ₀::Vector, T::Real, Nₜ::Int) = solveschrodinger(expmatrix(K, ψ₀, T/(Nₜ)), T, Nₜ-1)

function fftcₖ(ψ::Base.Callable, t::Real, Nₓ::Int, L::Real)
    X = range(-L/2,L/2,Nₓ)
    fft(apply(x->ψ(x,t), X))
end


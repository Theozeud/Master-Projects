using FFTW
using LinearAlgebra

# Base of L²
eₖ(k::Int, x::Real, L::Real) = exp(eltype(x)(2π) * im * k *x /L )

# Function to compute the Fourier coefficent of ψ : cₖ = <eₖ,ψ>
cₖ(ψ::Base.Callable, t::Real, k::Int, X::Discretization) = 1/length(X) * eₖ.(k,-X.points,X.sup-X.inf)'*ψ.(X.points,t)

# Function to compute the Fourier coefficent of any derivative of ψ : derivcₖ = <eₖ,∂ψ>
cₖ∂(order::Int, ψ::Base.Callable, t::Real, k::Int, X::Discretization) = (eltype(X)(2π)*im*k/(X.sup-X.inf))^order * cₖ(ψ, t, k, X)

approx(ψ::Base.Callable, X::Discretization, N::Int) = (x,t) -> mapreduce(k->eₖ(k,x,X.sup-X.inf)*cₖ(ψ, t, k, X),+, range(-N,N))

approx∂(order::Int, ψ::Base.Callable, X::Discretization, N::Int) = (x,t) -> mapreduce(k->eₖ(k,x,X.sup-X.inf)*cₖ∂(order,ψ, t, k, X),+, range(-N,N))


# Functions using FFT

function FFTapprox(ψ::Base.Callable, X::Discretization, t::Real)
    c̃ₖ = fft(ψ.(X,t))
    Ψ = ifft(c̃ₖ)
    Ψ,c̃ₖ
end

function FFTapprox∂(order::Int, ψ::Base.Callable, X::Discretization, t::Real)
    Nₓ = length(X)
    derivVector = (fftfreq(Nₓ,Nₓ)*2π/L*im).^order
    c̃ₖ = fft(ψ.(X.points,t))
    Ψ = ifft(derivVector .* c̃ₖ)
    Ψ,c̃ₖ
end

# Functions to solve the schrodinger equation

onestep(K::Matrix, ψ₀::Vector, t::Real) = exp(-im*K*t)*ψ₀
onestep(K::Matrix, V::Vector, ψ₀::Vector, t::Real) = exp(-im*K*t)*exp(-im*V*t)*ψ₀

function solveschrodinger(ψ₀::Vector, T::Real, Nₜ::Int, K::Matrix, V::Matrix)
    next_step = onestep(K, V[:,1], ψ₀, T/(Nₜ))
    [next_step, solveschrodinger(next_step, T, Nₜ-1, K, V[:,2:end])...]
end


function FFTsolveschrodinger(ψ₀::Vector, X::Discretisation, T::Real, Nₜ::Int, V)
    # Some parameters
    Nₓ = length(X)
    Δt = T/Nₜ

    # Potential
    V = hcat([V.(X,t) for i in range(Δt, T, Nₜ)]...)

    # Matrix to reach next each step
    stepΔ = exp(-im*eltype(ψ₀)(0.5)*fftfreq(Nₓ,Nₓ)*Δt)
    stepV = exp(-im*V*Δt)

    # Matrix to store the solution
    sol = zeros(Nₓ,Nₜ+1)
    sol[:,1]=ψ₀

    # Loop to compute the solution for each time step
    for i in 2:Nₜ+1
        sol[:,i] = ifft(stepΔ * fft(stepV[i] * sol[:,i-1]))
    end

    sol
end


function Kmatrix(order::Int, N::Int, X::Discretization)
    @show Nₓ = length(X)
    @show krange = collect(-N:N)
    @show power = krange*X.points'
    @show ematrix = Diagonal(krange.^(order/2)) * (fill(exp(2*im*π/(X.sup-X.inf)),2N+1, Nₓ)).^power
    @show factor = 1/Nₓ*(2*im*π/(X.sup-X.inf))^order*(-1)^(order/2)
    @show factor * ematrix' * ematrix[end:-1:1,:]
end




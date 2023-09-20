using FFTW
using Plots

apply(f::Base.Callable, X) = [f(x) for x in X]

eₖ(k::Int, L::Real, x::Real) = exp(2π * im * k *x /L )

function cₖ(ψ::Base.Callable, t::Real, Nₓ::Int, L::Real, k::Int)
    X = range(-L/2,L/2,Nₓ)
    1/Nₓ * mapreduce(n->eₖ(k,L,-X[n])*ψ(X[n],t),+, eachindex(X))
end

function computeApprox(ψ::Base.Callable, L::Real, Nₓ::Int, N::Int)
    (x,t) -> mapreduce(k->eₖ(k,L,x)*cₖ(ψ, t, Nₓ, L, k),+, range(-N,N))
end


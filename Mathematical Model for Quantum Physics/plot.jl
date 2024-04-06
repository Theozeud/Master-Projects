using Plots
using LaTeXStrings

function plotψ(ψ₁::Base.Callable, ψ₂::Base.Callable, X::Discretization, t::Real = 0)
    vX = collect(X.points)
    pltreal = plot(vX, real(ψ₁.(vX,t)), label = "Exact", size=(700,400), ylabel = L"Re(ψ)", line=(:black, 0.5, 6, :solid))
    plot!(pltreal, vX,real(ψ₂.(vX,t)), label = "Approximation", size=(700,400))
    
    pltimag = plot(vX, imag(ψ₁.(vX,t)), label = "Exact", size=(700,400), ylabel = L"Im(ψ)", line=(:black, 0.5, 6, :solid))
    plot!(pltimag, vX,imag(ψ₂.(vX,t)), label = "Approximation", size=(700,400))
    
    plot(pltreal, pltimag, layout=(1,2), plot_title="Approximation by the Fourier Transform")
    xlabel!(L"x")
    
end

function plotψ(ψ₁::Base.Callable, ψ₂X::Vector, _X, t::Real = 0)
    X = collect(_X.points)
    plot(X, real(ψ₁.(X,t)), label = "Exact")
    plot!(X,real.(ψ₂X), label = "Approximation")
end

function plotcₖ(C::Vector)
    Nₓ = length(C)
    plot(eachindex(C),real.(C))
    vline!([floor(Nₓ/2)])
end
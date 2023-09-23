using Plots

function plotψ(ψ₁::Base.Callable, ψ₂::Base.Callable, _X, t::Real = 0)
    X = collect(_X.points)
    plot(X, apply(x->real(ψ₁(x,t)), X), label = "Exact")
    plot!(X,apply(x->real(ψ₂(x,t)), X), label = "Approximation")
end

function plot_cₖ(C::Vector)
    plot(eachindex(C),C)
    vline!()

end
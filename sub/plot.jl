function plotψ(ψ₁::Base.Callable, ψ₂::Base.Callable, X, t::Real = 0)
    plot(X, apply(x->real(ψ₁(x,t)), X), label = "Exact")
    plot!(X,apply(x->real(ψ₂(x,t)), X), label = "Approximation")
end

function plot_cₖ()


end
using Distances
using Optim

function loss(K::Matrix, ψ::Vector, q::Vector, ε::Real)
    @assert length(collect(axes(K)[2])) == length(ψ)
    @assert length(collect(axes(K)[1])) == length(q)
    sqeuclidean(K*ψ,q)+ ε*norm(ψ.-.2)^2
end

defaultInit(ψ::Vector) = zero(ψ)


function optimizePSD(K::Matrix, q::Vector, ε::Real, method = LBFGS(), ψ₀ = defaultInit(length(q)), loss = loss)
    Optim.optimize(ψ->loss(K,ψ,q,ε), ψ₀;iterations = 100000)
end


function computeError(ψₙᵤₘ, ψₜᵣᵤₑ)
    sqeuclidean(ψₙᵤₘ,ψₜᵣᵤₑ)
end

function plotError(ψₙᵤₘ, ψₜᵣᵤₑ)
    
end

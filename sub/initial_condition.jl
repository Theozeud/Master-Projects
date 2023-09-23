using LinearAlgebra

# This files contains differents initials conditions

function gaussian(dis::Discretization, μ::Real = (dis.inf+dis.sup)/2, σ::Real = (dis.sup-dis.inf)/4)
    [1/sqrt(2π*σ)*exp(-(x-μ)^2/(2σ^2)) for x in dis]
end


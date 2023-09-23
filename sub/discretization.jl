struct Discretization{T}
    N::Int
    inf::Real
    sup::Real
    points::T
    function Discretization(N::Int,inf::Real,sup::Real, points)
        @assert N == length(points)
        @assert inf < sup
        new{typeof(points)}(N,inf,sup,points)
    end
end

uniformDiscretization(inf::Real,sup::Real, N::Int) = Discretization(N,inf,sup,collect(range(inf,sup,N+1))[1:N])

Base.iterate(d::Discretization, state::Int = 1) = state > d.N ? nothing : (d.points[state], state+1)
Base.eachindex(d::Discretization) = eachindex(d.points)
Base.getindex(d::Discretization, n::Int) = d.points[n]
Base.length(d::Discretization) = d.N
Base.eltype(d::Discretization) = eltype(d.points)
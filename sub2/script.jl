include("cld.jl")
include("plot.jl")
include("shape.jl")



shape = Cube(4)



#rot_shape = apply(p->Rotation(π,π/3,π/4, p),vertices(shape))

#plot2D(rot_shape)


@time cld = computeCLD(shape, 10000; Φ = 0, ϕ = 0)

plotCLD(cld)

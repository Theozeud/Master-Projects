include("../scr/shape.jl")
include("../scr/cld.jl")
include("../scr/plot.jl")

include("../scr/opti.jl")
include("../scr/random.jl")


# Inputs
sizeR = 10
sizeL = 1000
ntir = 1000000
ntirCLD = ntir

# Random size
random = Uniform(5, 10)
#random = Normal(5, 0.1)
#R = tir(random, sizeR)

# Given L (comes from the mesurement)
#L = range(1,4,sizeL) # Problème ??

# Base Shape
shape = Cube(1)

# L comes is randomly chosen in the experimental cld obtained with the shape
Ltir = computeCLD(shape, ntir)

# Utils function
#bins, bins_number = computeCumulCLD(shape, ntir, 10000)

# Computation of the matrix K
rMax = 20
lMax = maximum(Ltir)
rList = collect(range(0, rMax, sizeR))
lList = collect(range(0, lMax, sizeL))
cumulCLD = computeCumulCLD(Ltir)
K = matrixCLD(sizeL, sizeR, sorted(Ltir), rList, cumulCLD)

# Computation of the CLD for different size
function computeCLDdifferentSize()
    cld = zeros(ntir)
    # Random measurement
    for i in 1:ntir
        r = tir(random, 1)[1]
        cp = Cube(r)
        cld[i] = computeCL(cp)
    end
    bins, bins_number = computeCumulCLD(cld, 10000)
    [probaCLD(l, cumulCLD) for l in range(0, lMax, sizeL)]
end
q = computeCLDdifferentSize()

#function computeCLDdifferentSize(sizeR::Int, rList::Vector)
#    q = zeros(sizeR)
#    for r in rList
#        cld = 

# Optimisation
ψ₀ = [0.2 for i in 1:sizeR]
method = LBFGS()
ε = 0.1

opt = optimizePSD(K, q, ε, method, ψ₀)

#trueresult = [law(random, r) for r in R]

trueresult = apply(x -> x >= 5 && x <= 10 ? 0.2 : 0, rList)

res = Optim.minimizer(opt)

#plotOptiRes(rList, res)

qOpti = K * trueresult
plot(lList, qOpti)
plot!(lList, q)
using Plots
using LinearAlgebra
using LaTeXStrings
using FFTW
using Interpolations

function dynamics(ψ₀fun= x->exp(-x^2), Vfun = (x,t) -> 0; L=10, Nₓ=1000, T=4, Nₜ=1000)

    # Space discretization
    space = range(-L/2,L/2,Nₓ)[1:Nₓ]

    # Initial condition
    ψ₀ = ψ₀fun.(space)
    
    # Timestep
    Δt = eltype(ψ₀)(T/Nₜ)

    # Time discretization
    time = range(Δt,T,Nₜ)
    
    # Potential
    V = hcat([Vfun.(space,tᵢ) for tᵢ in time]...)

    # Matrix to reach next each step
    stepΔ = exp.(-im/2*eltype(ψ₀)(0.5)*(fftfreq(Nₓ,Nₓ)*2*π/L).^2*Δt)
    stepV = exp.(-im*V*Δt)

    # Matrix to store the solution
    sol = zeros(Complex, Nₓ,Nₜ+1)
    sol[:,1]=ψ₀
    
    # Loop to compute the solution for each time step
    for i in 2:Nₜ+1
        stepV[i] .* sol[:,i-1]
        fft(stepV[i] .* sol[:,i-1])
        stepΔ .* fft(stepV[i] .* sol[:,i-1])
        sol[:,i] = ifft(stepΔ .* fft(stepV[i] .* sol[:,i-1]))
    end
    sol
end

function animation(ψ::Matrix, duration::Real = 10, fps::Int = 30, L::Real = 10)

    @assert L > 0
    @assert duration > 0
    @assert fps > 0
    
    # set the min and maximum values of the plot, to scale the axis
    m = min(0, min(real.(ψ)...), min(imag.(ψ)...))
    M = max(abs.(ψ)...)

    # Total frame
    nbframe = fps * duration

    # space and time vector
    Nₓ,Nₜ = size(ψ)
    time = range(0,1,Nₜ)
    space = range(-L/2,L/2,Nₓ+1)[1:Nₓ]

    anim = @animate for n in 1:nbframe

        t = n / (fps*duration)

        linear_interpolation(time, ψ[1,:])
        linear_interpolation(time, ψ[1,:])(t)

        interp_linear = [linear_interpolation(time, ψ[i,:])(t) for i in eachindex(ψ[:,1])]
        
        pltreal = plot(space,real.(interp_linear), label = "Real", size=(700,400), xlim = [-L,L], ylim= [m,M])
        pltimag = plot(space,imag.(interp_linear), label = "Imag", size=(700,400), xlim = [-L,L], ylims= [m,M])
        pltabs  = plot(space, abs.(interp_linear), label = "Abs",  size=(700,400), xlim = [-L,L], y= [m,M])

        plot(pltreal, pltimag, pltabs, layout=(1,3), plot_title="Approximation by the Fourier Transform")

    end
    gif(anim, fps = fps)
end


@time animation(dynamics(Nₓ=100,T=10,Nₜ=100))
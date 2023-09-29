using Plots

#=
# Simulates the SchrÃ¶dinger dynamics iâˆ‚t = -1/2 Ïˆ'' + V(x,t) Ïˆ, with the pseudospectral method
# on an interval [-L,L] with periodic boundary conditions, with Nx grid points
# The simulation proceeds from 0 to T, with Nt time steps.
# Initial condition is given by the function psi0_fun(x), potential by the function V_fun(x,t)
# Returns an array psi[ix, it]
def dynamics(psi0_fun=(lambda x: np.exp(-x**2)), V_fun=(lambda x,t: 0), L=10, Nx=1000, T=4, Nt=1000):
    # TO DO


# Plots the return value psi of the function "dynamics", using linear interpolation
# The whole of psi is plotted, in an animation lasting "duration" seconds (duration is unconnected to T)
# L argument is only for x axis labelling
def plot_psi(psi, duration=10, frames_per_second=30, L=10):
    
    

    # set the axis once and for all
    ax.set(xlim=[-L,L], ylim=[m,M], xlabel='x', ylabel='psi')
    # dummy plots, to update during the animation
    real_plot = ax.plot(x_data, np.real(psi[:, 0]), label='Real')[0]
    imag_plot = ax.plot(x_data, np.imag(psi[:, 0]), label='Imag')[0]
    abs_plot  = ax.plot(x_data, np.abs(psi[:, 0]), label='Abs')[0]
    ax.legend()

    # define update function as an internal function (that can access the variables defined before)
    # will be called with frame=0...(duration*frames_per_second)-1
    def update(frame):
        print(frame)
        # get the data by linear interpolation
        t = frame / (duration * frames_per_second)
        psi_t = np.array([np.interp(t, t_data, psi[i, :]) for i in range(np.size(psi,0))])
        # update the plots
        real_plot.set_ydata(np.real(psi_t))
        imag_plot.set_ydata(np.imag(psi_t))
        abs_plot.set_ydata(np.abs(psi_t))

    ani = animation.FuncAnimation(fig=fig, func=update, frames=duration*frames_per_second, interval=1000/frames_per_second)
    return ani

=#

include("discretization.jl")

function animation(ψ, Nₜ, X::Discretization, T::Real, fps::Int = 30)
    

    # set the min and maximum values of the plot, to scale the axis
    #m = min(0, min(real.(ψX)...), min(imag.(ψX)...))
    #M = max(abs.(ψX)...)
    
    time = range(0, T, Nₜ)
    vX = collect(X.points)

    anim = @animate for i in 1:Nₜ
        
        pltreal = plot(vX,real.(ψ.(vX,time[i])), label = "Real", size=(700,400))
        pltimag = plot(vX,imag.(ψ.(vX,time[i])), label = "Imag", size=(700,400))
        pltabs  = plot(vX, abs.(ψ.(vX,time[i])), label = "Abs",  size=(700,400))

        plot(pltreal, pltimag, pltabs, layout=(1,3), plot_title="Approximation by the Fourier Transform")

    end
    gif(anim, fps = fps)
end
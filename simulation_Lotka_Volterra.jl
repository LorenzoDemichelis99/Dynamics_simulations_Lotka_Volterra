# this file implements the simulation of the dynamics

# The function implementing the simulation takes the following parameters:
# G: the underlying graph of the dynamics
# T: the number of time steps
# Δ: the time step

function simulation_Lotka_Volterra(G::SimpleGraph, T::Int64, Δ::Float64, x_0::Function, m::Float64, γ::Float64, σ::Float64, K::Int64)
    # initialization
    sim = simulation_data_Lotka_Volterra(G, T, Δ, m, γ, σ, K, x_0)
    N = nv(G)
    logx = 0.0

    # setting up the progress meter
    p = Progress(T; showspeed=true)

    # simulation of the dynamics
    @inbounds for t in 2:T
        for i in 1:N
            logx = log(sim.x[i,t-1]) + (Δ * (1 - sim.x[i,t-1]))
            for j in G.fadjlist[i]
                logx += (Δ * sim.α[i,j] * sim.x[j,t-1])
            end

            sim.x[i,t] = exp(logx)

            logx = 0.0
        end
        next!(p)
    end

    return sim
end
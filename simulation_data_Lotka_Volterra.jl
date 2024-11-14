# this file is used to define the struct which contains all the important information about the simulation of the dynamics of the s-spin model

# The important quantities for the simulations are:
# G: underlying graph of the system
# α: matrix of the coupling constants
# T: number of time steps of the simulation
# Δ: time steps
# x: matrix containing the trajectories generated in the simulation
# m: parameter m of the distribution of disorder
# γ: parameter γ of the distribution of disorder
# σ: parameter σ of the distribution of disorder
# K: average connectivity of the graph

mutable struct simulation_data_Lotka_Volterra
    G::SimpleGraph
    α::Matrix{Float64}
    T::Int
    Δ::Float64
    x::Matrix{Float64}
    m::Float64
    γ::Float64
    σ::Float64
    K::Int64
end

# this is the constructor of the struct defined above
function simulation_data_Lotka_Volterra(G::SimpleGraph, T::Int64, Δ::Float64, m::Float64, γ::Float64, σ::Float64, K::Int64, x0_init::Function)
    # extraction of the number of vertices
    n = nv(G)

    # allocation of the required memory
    α = Matrix{Float64}(undef, (n,n))
    x = Matrix{Float64}(undef, (n,T))
    
    # creation of the distribution of the couplings
    μ = [0, 0]
    Σ = [1.0 γ; γ 1.0]
    z = MvNormal(μ, Σ)

    # initialization of the required quantities
    coupling!(G, α, m, K, σ, z)
    init_x!(x, x0_init)

    # initialization of the struct
    simulation_data_Lotka_Volterra(G, α, T, Δ, x, m, γ, σ, K)
end

# this is the function that builds the matrix of couplings
function coupling!(G::SimpleGraph, α::Matrix{Float64}, m::Float64, K::Int64, σ::Float64, z::H) where H <: AbstractMvNormal
    α .= 0.0
    α_1 = 0.0
    α_2 = 0.0
    @inbounds for i in 1:size(α, 1)
        α[i,i] = 0.0
        for j in i+1:size(α, 2)
            if has_edge(G,i,j)
                α_1, α_2 = rand(z) 
                α[i,j] = (m/K) + ((σ/sqrt(K))*α_1) 
                α[j,i] = (m/K) + ((σ/sqrt(K))*α_2) 
            else
                α[i,j] = 0.0
            end
        end
    end
end

# this is the function that initializes the trajectories of the degrees of freedom
function init_x!(x::Matrix{Float64}, x0_init::Function)
    x .= 0.0
    @inbounds @Threads.threads for i in 1:size(x,1)
        x[i,1] = x0_init(x[i,1])
        x[i,2:end] .= 0.0
    end
end
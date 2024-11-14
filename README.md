# Simulation of the dynamics of the random Lotka-Volterra model on a generic graph

## Structure of the files

This repository contains a Jupyter notebook and some Julia files to implement the simulation of the random Lotka-Volterra model on a generic graph. 
The file structure is the following: simulation_data_Lotka_Volterra.jl contains the code that defines the mutable struct used to store all the important information about the simulation and the code, together with the constructor of the struct and the functions that initialize the relevant quantities of the simulation; simulation_Lotka_Volterra.jl contains the function implementing the actual simulation; simulation_tools_Lotka_Volterra.jl contains all the utility functions required to plot and save the results of the simulation.


## Dynamics implemented on the graph

Let $G=(V,E)$ be the graph on which the dynamics is defined, with $V$ being the node set and $E$ being the edge set; the evolution of the degree of freedom $x_{i}(t)$, with $i \in V$, is determined by the following stochastic differential equation:
```math
\begin{equation}
    \frac{dx_{i}(t)}{dt} = x_{i}(t) \bigg( 1 - x_{i}(t) + \sum_{j \in \partial i} \alpha_{ij} x_{j}(t) \bigg)
\end{equation}
```
with $\alpha_{ij}$ being the coupling constant between node $i$ and node $j$.

For what concerns the coupling constants, it is possible to introduce an asymmetry in the interactions by considering the pair $(\alpha_{ij}, \alpha_{ji})$ as independent random variables such that:
```math
\begin{equation}
    \alpha_{ij} = \frac{m}{K} + \frac{\sigma}{\sqrt{K}} z_{ij}
\end{equation}
```
where $\langle z_{ij} \rangle = 0$, $\langle z_{ij}^{2} \rangle = \langle z_{ji}^{2} \rangle = 1$ and $\langle z_{ij} z_{ji} \rangle = \gamma$.

The equations describing the evolution of the degrees of freedom are integrated using the Euler-Maruyama scheme, where time is discretized as $t=n\Delta$, with $\Delta$ being a small time step. By setting $x_{i}^{n} = x(t=n\Delta)$ the equations of the dynamics in discrete time and in log-space becomes:
```math
\begin{equation}
    \log x_{i}^{n+1} = \log x_{i}^{n} + \Delta \bigg( 1 - x_{i}^{n} + \sum_{j \in \partial i} \alpha_{ij} x_{j}^{n}\bigg)
\end{equation}
```

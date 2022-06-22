"""
One photon moving in Minkowski spacetime
"""

include("metrics/minkowski.jl")
using DifferentialEquations
using Plots
using .Minkowski

# Initial Condition
x0 = [5., pi/2 , 0., 0.]        # initial position
k0 = [0., -0.1, 1., 0.]         # initial momentum (except kt = k0[4])
k0[4] = sqrt(k0[1]^2 + x0[1]^2 * k0[2]^2 + x0[1]^2*sin(x0[2])^2 *k0[3]^2)
# normalization condition for the momentum

# Lower the indices of the momentum
k_0 = zeros(4)
g = metric(x0)
k_0[1] = g[1]*k0[1]
k_0[2] = g[2]*k0[2]
k_0[3] = g[3]*k0[3] + g[5]*k0[4]
k_0[4] = g[4]*k0[4] + g[5]*k0[3] 

# Complete initial state
xk0 = [x0; k_0]


# Solving the geodesic equations 
prob = ODEProblem(geodesics, xk0, (0., 5.), saveat=(1:1:5))
sol = solve(prob)


# Introduce cartesian coordinates
x = sol[1,:].* sin.(sol[2,:]) .* cos.(sol[3,:])
y = sol[1,:].* sin.(sol[2,:]) .* sin.(sol[3,:])
z = sol[1,:].* cos.(sol[2,:])


# Plot the trajectory

pl = plot(x,y,z, #seriestype=:scatter, 
    xlims= (-10,10),
    xlabel="x", ylabel="y", zlabel="z")

savefig(pl, "onephoton.pdf")
display(pl)

readline()


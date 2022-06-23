include("screen/imageplane.jl")
include("metrics/minkowski.jl")

using .Minkowski
using .ImagePlane

using DifferentialEquations
using Plots
using Printf



"""
geoInteg(geodesics, xk0)
------------------------------------------------------------
Integrates the geodesics until the photon reaches the 
equatorial plane
------------------------------------------------------------
Arguments:
geodesics : function containing the geodesic equations
xk0 : Vector with the initial values of position and 
      momentum at the image plane.
      [r, theta, phi, t, k_r, k_theta, k_phi, k_t]
------------------------------------------------------------
Returns:
xk0 (modified) : Vector with the final state of the photon
			[r, theta, phi, t, k_r, k_theta, k_phi, k_t]
------------------------------------------------------------
"""
function geoInteg(geodesics, xk0)
	# step of integration (negative to go backwards in time)
    dlmbda = -.005
    # counter of integration steps
    j = 1
    # Max value of the counter to stop the process
    jmax = 1000000

    # Main loop verifying the position of the photon
    while xk0[1]*cos(xk0[2]) > 1e-5 && j<jmax
        prob = ODEProblem(geodesics, xk0, (0., dlmbda))
        sol = solve(prob, save_everystep=false)
        xk0 = sol[1:8,end]
        j +=1
    end
    # Returns the final state or empty if the photon does
    # not reach the equatorial plane
    if j<jmax
        return xk0
    else
        return[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]
    end
end


# Initial parameters of the system 
# (These will be loaded from a configuration file in a 
#  future version)
screensize = 20
N = 50
D = 100.
inclntn = pi/4


# Screen definition
alphaRange, betaRange, numPix = screen(screensize, N)
totalPix = numPix*numPix

# Final radial position of the photon
rfdata = zeros(length(alphaRange), length(betaRange))

# Integration for all the photons at the image plane
k=0
for i in 1:length(alphaRange), j in 1:length(betaRange)
    xkf = geoInteg(geodesics, initConds(alphaRange[i],betaRange[j], D, inclntn, metric))
    rfdata[i,j] = xkf[1]
    @printf("%.2f %%", k*100/totalPix)
    println()
    global k += 1
end


# Image
img = heatmap(1:size(rfdata,1), 1:size(rfdata,2), rfdata)
savefig(img, "diskimage.pdf")
display(img)
readline()
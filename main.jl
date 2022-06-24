
include("config.jl")
import .config as cfg

println("\nReading Parameters \n")


if cfg.st == "Minkowski"
	include("metrics/minkowski.jl")
	using .Minkowski
	println("\nSpacetime: Minkowski")
else
	println("\nThe selected spacetime is not available\n")
	exit()
end

if cfg.structure == 1
	include("structures/simpledisk.jl")
	using .SimpleDisk
	r_in = cfg.r_in
	r_out = cfg.r_out
	println("Accretion Structure: Simple Disk")
else
	println("\nThe selected accretion structure is not available\n")
	exit()
end

if cfg.screenType == 1
	include("screen/imageplane.jl")
	using .ImagePlane
	distance = cfg.distance
	inclination = cfg.inclination
	screensize = cfg.screensize
	numberOfPixels = cfg.numberOfPixels
	println("Screen Type: Image Plane")
else
	println("\nThe selected screen is not available\n")
	exit()
end


using DifferentialEquations
using Printf
#using BenchmarkTools
using JLD



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



# Screen parameters definition
alphaRange, betaRange, numPix = screen(screensize, numberOfPixels)
totalPix = numPix*numPix


# Energy information
photonEnergy = zeros(numPix, numPix)


k=0 # Percentage of advance counter
# Integration for all the photons at the image plane
@time for i in 1:length(alphaRange), j in 1:length(betaRange)
    xkf = geoInteg(geodesics, initConds(alphaRange[i],betaRange[j], distance, inclination, metric))
    photonEnergy[i,j] = emission(xkf[1]; r_in=r_in, r_out=r_out)
    @printf("%.2f %%", k*100/totalPix)
    println()
    global k += 1
end

# Save the image data
save("data/ImageData.jld", "ImageData", photonEnergy)



if cfg.showImage == true 
	include("common/image.jl")
	import .Image as img
	img.createImage(photonEnergy, numPix, cfg.saveImage; cfg.imageName)
else
	println("\nImage will not be saved.")
end



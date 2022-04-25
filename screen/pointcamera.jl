module PointCamera

export screen, Photon, initCoords

"""
screen(FoV, N::Integer)
------------------------------------------------------------
Defines a square NxN pixels screen.
Returns the ranges of the cordinates alpha and beta 
in the image plane.
------------------------------------------------------------
Arguments:
FoV : Field of View in radians
N : Number of pixels wanted
------------------------------------------------------------
"""
function screen(FoV, N::Integer)
	# We will use an odd number of pixels in each direction
	if iseven(N) 
		numPixels = N + 1
	else
		numPixels = N
	end

	alphaRange = LinRange(-FoV/2, FoV/2, numPixels)
    betaRange = LinRange(-FoV/2, FoV/2, numPixels)
    println("\nPoint Camera")
    println("\nSize of the screen in Pixels: $numPixels X $numPixels")
    println("\nTotal Number of Pixels: $(numPixels*numPixels)\n")
    return alphaRange, betaRange, numPixels
end


"""
Photon
------------------------------------------------------------
Structure to defibe the set of photons
------------------------------------------------------------
Attributes:
alpha, beta : Coordinates of the photon in the image plane
rf : Will store the location of the photon at the equatorial 
     plane of the accretion structure 
------------------------------------------------------------
"""
mutable struct Photon
    alpha::Float64
    beta::Float64
    rf::Float64
end


"""
initCoords(p::Photon, D, inclntn, metric; K0 = 1.)
------------------------------------------------------------
Calculates the initial parameters for a photon
------------------------------------------------------------
Arguments:
p : Photon 
D : Distance between compact object and point camera [kpc]
inclntn : Inclination angle [radians]
metric : metric tensor in the region of the camera. ::Vector
		 grr = g[1]
		 gthth = g[2]
		 gphph = g[3]
		 gtt = g[4]
		 gtphi = g[5]
K0 : Magnitude of the momentum vector for the photon
------------------------------------------------------------
Returns:
xk0 : Vector with the initial values of position and 
      momentum at the camera location.
      [r, theta, phi, t, k_r, k_theta, k_phi, k_t]
------------------------------------------------------------
"""
function initCoords(p::Photon, D, inclntn, metric; K0 = 1.)
	# The camera is a point, therefore the initial position is
	# the same for all photons!
	xk0 = zeros(8)
	xk0[1] = D
    xk0[2] = inclntn
    xk0[3] = 0.
    xk0[4] = 0.
        
    # Initial 4-momentum components
    kr =  K0*cos(p.beta)*cos(p.alpha)
    ktheta = -(K0*sin(p.beta)/D)
    kphi = K0*cos(p.beta)*sin(p.alpha)/(D*sin(inclntn))
    kt = sqrt(kr^2 + xk0[1]^2 * ktheta^2 + xk0[1]^2*sin(xk0[2])^2 *kphi^2)
    
    g = metric(xk0[1:4])
    xk0[5] = g[1]*kr
    xk0[6] = g[2]*ktheta
    xk0[7] = g[3]*kphi + g[5]*kt
    xk0[8] = g[4]*kt + g[5]*kphi

    return xk0
end

end # end of module


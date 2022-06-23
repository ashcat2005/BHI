module ImagePlane

export screen, Photon, initConds

"""
------------------------------------------------------------
This module defines the image plane as screen to receive the 
photons and the initial conditions for each photon
------------------------------------------------------------
Functions:
screen: defines a square NxN pixels screen
Photon: structure to create each of the photons that will be 
        traced back onto the accretion structure
initCoords: defines the initial coordinates for a Photon
------------------------------------------------------------
@author: ashcat - 2022
"""


"""
screen(Ssize, N::Integer)
------------------------------------------------------------
Defines a square NxN pixels screen.
------------------------------------------------------------
Arguments:
Ssize: Half of the Screen size (in units of M). The screen 
       will go from -Ssize to +Ssize
N: number of pixels in the screen
------------------------------------------------------------
Returns: 
The ranges of the celestial cordinates and the total number 
of pixels in the screen
Alpha : Horizontal celestial coordinate
Beta : Vertical celestial coordinate
numPixels : Number of pixels (odd number) in each direction
------------------------------------------------------------
"""
function screen(Ssize, N::Integer)
	# We will use an odd number of pixels in each direction
	if iseven(N) 
		numPixels = N + 1
	else
		numPixels = N
	end

	alphaRange = LinRange(-Ssize, Ssize, numPixels)
    betaRange = LinRange(-Ssize, Ssize, numPixels)
    println("\nImage Plane")
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
FIRST FORM
initConds(p::Photon, D, inclntn, metric; K0 = 1.)
------------------------------------------------------------
Calculates the initial parameters for a photon
------------------------------------------------------------
Arguments:
p : Photon 'structure'
D : Distance between compact object and image plane [kpc]
inclntn : Inclination angle [radians]
metric : metric tensor in the region of the camera. ::Vector
		 grr = g[1]
		 gthth = g[2]
		 gphph = g[3]
		 gtt = g[4]
		 gtphi = g[5]
K0 = 1 : Magnitude of the momentum vector for the photon
------------------------------------------------------------
Returns:
xk0 : Vector with the initial values of position and 
      momentum at the image plane.
      [r, theta, phi, t, k_r, k_theta, k_phi, k_t]
------------------------------------------------------------
"""
function initConds(p::Photon, D, inclntn, metric; K0 = 1.)
	# Transformation from (Alpha, Beta, D) to (r, theta, phi)
	xk0 = zeros(8)
	xk0[1] = sqrt(p.alpha^2 + p.beta^2 + D^2)
    xk0[2] = acos((p.beta*sin(inclntn) + D*cos(inclntn))/xk0[1])
    xk0[3] = atan(p.alpha/(D*sin(inclntn) - p.beta*cos(inclntn)))
    xk0[4] = 0.
        
    # Initial 4-momentum components
    kr =  (D/xk0[1])*K0

    aux = p.alpha^2 + (-p.beta*cos(inclntn) + D*sin(inclntn))^2

    ktheta = (K0/sqrt(aux))*(-cos(inclntn) 
                           +(p.beta*sin(inclntn) + D*cos(inclntn))
                *(D/(xk0[1]^2)))

    kphi = -p.alpha*sin(inclntn)*K0/aux
    kt = sqrt(kr^2 + xk0[1]^2 * ktheta^2 + xk0[1]^2*sin(xk0[2])^2 *kphi^2)
    
    g = metric(xk0[1:4])
    xk0[5] = g[1]*kr
    xk0[6] = g[2]*ktheta
    xk0[7] = g[3]*kphi + g[5]*kt
    xk0[8] = g[4]*kt + g[5]*kphi

    return xk0
end


"""
SECOND FORM
initConds(alpha, beta, D, inclntn, metric; K0 = 1.)
------------------------------------------------------------
Calculates the initial parameters for a photon
------------------------------------------------------------
Arguments:
alpha : Horizontal celestial coordinate of the photon
beta  : Vertical celestial coordinate of the photon
D : Distance between compact object and image plane [kpc]
inclntn : Inclination angle [radians]
metric : metric tensor in the region of the camera. ::Vector
         grr = g[1]
         gthth = g[2]
         gphph = g[3]
         gtt = g[4]
         gtphi = g[5]
K0 = 1 : Magnitude of the momentum vector for the photon
------------------------------------------------------------
Returns:
xk0 : Vector with the initial values of position and 
      momentum at the image plane.
      [r, theta, phi, t, k_r, k_theta, k_phi, k_t]
------------------------------------------------------------
"""
function initConds(alpha, beta, D, inclntn, metric; K0 = 1.)
    # Transformation from (alpha, beta, D) to (r, theta, phi)
    xk0 = zeros(8)
    xk0[1] = sqrt(alpha^2 + beta^2 + D^2)
    xk0[2] = acos((beta*sin(inclntn) + D*cos(inclntn))/xk0[1])
    xk0[3] = atan(alpha/(D*sin(inclntn) - beta*cos(inclntn)))
    xk0[4] = 0.
        
    # Initial 4-momentum components
    kr =  (D/xk0[1])*K0

    aux = alpha^2 + (-beta*cos(inclntn) + D*sin(inclntn))^2

    ktheta = (K0/sqrt(aux))*(-cos(inclntn) 
                           +(beta*sin(inclntn) + D*cos(inclntn))
                *(D/(xk0[1]^2)))

    kphi = -alpha*sin(inclntn)*K0/aux
    kt = sqrt(kr^2 + xk0[1]^2 * ktheta^2 + xk0[1]^2*sin(xk0[2])^2 *kphi^2)
    
    g = metric(xk0[1:4])
    xk0[5] = g[1]*kr
    xk0[6] = g[2]*ktheta
    xk0[7] = g[3]*kphi + g[5]*kt
    xk0[8] = g[4]*kt + g[5]*kphi

    return xk0
end


end # end of module


module Minkowski
export metric, geodesics, rEH
"""
------------------------------------------------------------
Module with the information for the Minkowskian metric in
spherical coordinates

ds^2 = -dt^2 + dr^2 + r^2 dtheta^2 + r^2 sin^2 (theta) dphi^2
------------------------------------------------------------
Properties:
- No curvature
- No event horizon
- No ISCO
------------------------------------------------------------
@author: ashcat - 2022
"""



"""
metric(x::Vector)
------------------------------------------------------------
Returns the metric components for the Minkowski metric
------------------------------------------------------------
Arguments:
The components of x are interpreted as
r = x[1]
theta = x[2]
phi = x[3]
t = x[4]
------------------------------------------------------------
Returns:
The metric components returned are
g[1] = grr
g[2] = gthth
g[3] = gphph
g[4] = gtt
g[5] = gtphi
------------------------------------------------------------
"""
function metric(x::Vector)
	grr = 1.  					
	gthth = x[1]^2				
	gphph = x[1]^2*sin(x[2])^2	
	gtt = -1.					
	gtphi = 0.  				
	return [grr, gthth, gphph, gtt, gtphi]
end



"""
geodesics(xk::Vector,p,t)
------------------------------------------------------------
Arguments:
Contains the geodesic equations in Hamiltonian form for the 
Minkowski metric. The components of xk are interpreted as
r = xk[1]
theta = xk[2]
phi = xk[3]
t = xk[4]
k_r = xk[5]
k_th = xk[6]
k_phi = xk[7] = L
k_t = xk[8] = -E

** The arguments p and t are not used (but needed by the 
integration algorithm)
------------------------------------------------------------
Returns:
Vector with the RHS of the geodesic equations
------------------------------------------------------------
"""
function geodesics(xk::Vector, p, t)
	sinth = sin(xk[2])
	drdl = xk[5]  						# dr/dlmbda = k_r
	dthdl = xk[6]/xk[1]^2				# dth/dlmbda = k_th/r^2
	dphdl = xk[7]/(xk[1]^2*sinth^2)		# dph/dlmbda = L/r^2sin(th)^2
	dtdl = -xk[8]						# dt/dlmbda = E

	dkrdl = xk[6]^2/xk[1]^3 + xk[7]^2/(xk[1]^3*sinth^2)
	# dkr/dlmbda = kth^2/r^3 + L^2/r^3*sin(th)^2
	dkthdl = (cos(xk[2])/sinth^3)*(xk[7]/xk[1])^2
	# dkth/dlmbda = cos(th)/sin(th)^3 * L^2/r^2
	dkphdl = 0.							# dkph/dlmbda = 0
	dktdl = 0.							# dkt/dlmbda = 0
	return [drdl, dthdl, dphdl, dtdl, dkrdl, dkthdl, dkphdl, dktdl] 
end



"""
rEH()
------------------------------------------------------------
Arguments:
None
------------------------------------------------------------
Returns:
Minkowski's metric has no horizon. Hence it returns the 
value 0.
------------------------------------------------------------
"""
function rEH()
	return 0.
end




end # end of module
module Minkowski
export metric, geodesics

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
	g = zeros(5)
	g[1] = 1.  					# grr
	g[2] = x[1]^2				# gthth
	g[3] = x[1]^2*sin(x[2])^2	# gphph
	g[4] = -1.					# gtt
	g[5] = 0.  					# gtphi
	return g
end



"""
geodesics(xk::Vector,p,t)
------------------------------------------------------------
Arguments
Contains the geodesic equations in Hamiltonian form for the 
Minkowski metric. The components of xk are
interpreted as
r = xk[1]
theta = xk[2]
phi = xk[3]
t = xk[4]
k_r = xk[5]
k_th = xk[6]
k_phi = xk[7] = L
k_t = xk[8] = -E

The arguments p and t are not used (but needed by the 
integration algorithm)
------------------------------------------------------------
Returns:
Vector with the RHS of the geodesics
------------------------------------------------------------
"""
function geodesics(xk::Vector, p, t)
	f = zeros(8)
	f[1] = xk[5]  						# dr/dlmbda = k_r
	f[2] = xk[6]/xk[1]^2				# dth/dlmbda = k_th/r^2
	f[3] = xk[7]/(xk[1]^2*sin(xk[2])^2)	# dph/dlmbda = L/r^2sin(th)^2
	f[4] = -xk[8]						# dt/dlmbda = E

	# dkr/dlmbda = kth^2/r^3 + L^2/r^2sin(th)^2
	f[5] = xk[6]^2/xk[1]^3 + xk[7]^2/(xk[1]^3*sin(xk[2])^2)
	# dkth/dlmbda = cos/th)/sin(th)^3 * L^2/r^2
	f[6] = (cos(xk[2])/sin(xk[2])^3)*(xk[7]/xk[1])^2
	f[7] = 0.							# dkph/dlmbda = 0
	f[8] = 0.							# dkt/dlmbda = 0
	return f 
end

end # end of module
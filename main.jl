

using DifferentialEquations
#using Plots

include("metrics/minkowski.jl")
include("screen/pointcamera.jl")

using .Minkowski
using .PointCamera


function geoInteg(geodesics, xk0)
	dlmbda = -.005
	j = 1
	jmax = 1000000
	while xk0[1]*cos(xk0[2]) > 1e-5 && j<jmax
		prob = ODEProblem(geodesics, xk0, (0., dlmbda))
		sol = solve(prob, save_everystep=false)
		xk0 = sol[1:8,end]
		j +=1
	end
	if j<jmax
		return xk0
	else
		return[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]
	end
end


N = 10
FoV = pi/4

D = 100.
inclntn = pi/4

alphaRange, betaRange, numPixels = screen(FoV, N)

photons = [Photon(a,b, 0.) for a in alphaRange, b in betaRange ] 

@time begin
#Threads.@threads for p in photons[1:20]
for p in photons
	xkf = geoInteg(geodesics, initCoords(p, D, inclntn, metric))
	p.rf = xkf[1]
	println(p.rf)
end
end





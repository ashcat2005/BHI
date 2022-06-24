module Verlet

export velocityVerlet

function velocityVerlet(ODE::Function, q0::Vector, lmbda::Vector)
	N = length(lmbda) # Independent parameter steps
	dlmbda = lmbda[2] - lmbda[1]
	q = zeros(N, 6)
	q[1,:] = q0
	for i in 1:N-1 
		v_half = q[i,4:6] + 0.5*ODE(q[i,:])[4:6]*dlmbda
		q[i+1,1:3] = q[i,1:3] + v_half*dlmbda
		q[i+1,4:6] = v_half + 0.5*ODE(q[i,:])[4:6]*dlmbda
	end
	return q
end

end # end of module
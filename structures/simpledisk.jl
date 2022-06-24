module SimpleDisk
export emission

"""
------------------------------------------------------------
This module defines a simple thin disk structure in the 
equatorial plane
------------------------------------------------------------
Functions:
structure: defines a simple disk structure in the equatorial 
           plane with a non-physical energy spectrum
------------------------------------------------------------
@author: ashcat - 2022
"""


"""
emission(r; r_in = 6., r_out = 10.)
------------------------------------------------------------
Defines a simple disk structure  with a non-physical 
energy spectrum proportional to the inverse of the radial
coordinate.
------------------------------------------------------------
Arguments:
r: radial coordinate in the disk surface
r_in : innner radius of the disk. Default = 6.
r_out : outer radius of the disk. Default = 10.
------------------------------------------------------------
Returns: 
Local flux of radiative energy produced at r
------------------------------------------------------------
"""
function emission(r; r_in = 6., r_out = 10.)
    if r >= r_in && r<= r_out
        return 1/r
    else
        return 0.
    end
end


end # end of module


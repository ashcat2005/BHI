module config


"""
------------------------------------------------------------
Configuration File    BHI
------------------------------------------------------------
@author: ashcat - 2022
"""



"""
Spacetime Selection
------------------------------------------------------------
Options and Parameters
- "Minkowski" º
   None
- "Schwarzschild"
   M : Mass of the black hole
- "Kerr"
   M : Mass of the black hole
   a : Spin parameter of the black hole
------------------------------------------------------------
"""
############################################################
st = "Minkowski"
M = 1.  # Default value just to define the length units 
        # (Not an actual mass!)
############################################################
# st = "Schwarzschild"
# M = 1.
############################################################
# st = "Kerr"
# M = 1.
# a = 0.5
############################################################




"""
Accretion Structure 
------------------------------------------------------------
Options and Parameters:

# 1: Simple Accretion Disk
     r_in
     r_out

# 2: Infinite Accretion Disk with a decreasing exponential 
     spectrum
     r_in
     r_out
     tau
     corotation

# 3: Novikov-Thorne Thin Accretion Disk (not yet!)
     r_in
     r_out
     corotation
------------------------------------------------------------
"""
############################################################
structure = 1
r_in = 6*M
r_out = 12*M
# tau
# corotation = true
############################################################



"""
Screen Definition
------------------------------------------------------------
Options and Parameters:

# 1: Image Plane
     distance
     inclination
     screensize
     numberOfPixels

# 2: Point Camera 
     Distance
     FieldOfView
------------------------------------------------------------
"""
############################################################
screenType = 1
distance = 100 # in [kpc]
inclination = pi/4
screensize = 10
numberOfPixels = 5
############################################################



"""
Image Visualization
------------------------------------------------------------
Options and Parameters:

showImage

saveImage
------------------------------------------------------------
"""
############################################################
showImage = true 
saveImage = true
imageName = "data/diskImage"

end # end of module

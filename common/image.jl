module Image 
export createImage


using Plots
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
function createImage(photonEnergy, numPix, saveImage; imageName="diskimage")
	img = heatmap(1:numPix, 1:numPix, photonEnergy)
	if saveImage==true
		savefig(img, imageName*string(numPix*numPix)*"px.png")
	end

	display(img)
	readline()
end

end # end of module

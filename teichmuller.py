from time import time
#from mpmath import mp  # Use the mpmath module to work to arbitary precision.

from surfaces import PantsDecomposition
from tiling import Tiling

# Start by building a pants decomposition of a surface.
# Here we are building one from 2 pairs of pants as 2 sets of
# cuff lengths are given. These are numbered pants 0 and pants 1 and each have
# 3 cuffs, numbered 0, 1 and 2. Each cuff length is given as a tuple and, in
# this case, they all have the same length, namely 1.0.
# Guing information has also been specified. These are tuples of the form:
#   (source_pant_index, source_pant_cuff, target_pant_index, target_pant_cuff, torsion).
# So:
#   (0,1,0,2,0.0)
# Should be read as:
#   Glue pants 0's cuff 1 to pants 0's cuff 2 with a 0.0 twist.
# And similarly:
#   (0,0,1,0,0.1)
# Should be read as:
#   Glue pants 0's cuff 0 to pants 1's cuff 1 with a 0.1 twist.
# Doing this now is optional and gluings can always be specified later using the
# glue_cuffs and unglue_cuffs functions of S.

#                          V The cuff lengths.
S = PantsDecomposition( ((1.0,4.0,1.0), (1.0,1.0,1.0)), \
                         ((0,1,0,2,0.0), (0,0,1,0,0.1), (1,1,1,2,-0.3)) )
#                          ^ The gluing information.

# So in this case the these gluings create a genus 2 surface with no boundary. You should check
# this by drawing a picture.

# We can double check some information about the surface we created just to make sure it
# is the right one.
genus, boundary_components = S.information()
print('This is a surface of genus %d with %d boundary components.' % (genus, boundary_components))
# This should be genus 2 with no boundary components.

# We now construct a 1000 right angled hexagon tile tiling of H^2 using S.
T = Tiling(S, number=1000)
# Note in general the tilings out to a certain radius look 'nicer'. So we could have done:
# T = Tiling(S, radius=7.0)

# If we want to look at the actual tiles involved we can either get them from the suface,
H = S.hexagon_decomposition()
# Or ask the tiling for the ones it used.
H = T.hexagons
# However, there should be very few reasons that you would ever want to do this.

# We can now ask for properties of this tiling, such as its length spectrum.
L = T.length_spectrum()

# The systole is the curve on S of shortest length. So we can determine its length by asking
# for the minimum length in the spectrum.
print 'The systole of this surface has length %f.' % min(L)

# The spectrum is returned in a random order, so a sensible thing to do is to first sort it.
L.sort()

# We could now look at some of the shorter curves.
print('The shortest 10 curves in its length spectrum have lengths:')
for l in L[:10]: print l

# But if we are only interested in seening if the spectrum contains a curve of length
# at most 2.2 for example, we can use the threshold function to stop contructing the
# spectrum early.
L1 = T.length_spectrum(threshold=2.2)

# Then there is a curve of length at most 2.2 iff the last element of the L1 list is
# at most 2.2. So we can do:
if L1[-1] < 2.2:
    print('There is a curve in the length spectrum with length < 2.2.')
else:
    print('There is NOT a curve in the length spectrum with length < 2.2.')

# But of course we could just use the built in function:
if T.epsilon_short_curve(2.2):
    print('Yes, there really is a curve with length < 2.2.')
else:
    print('No, there really is not a curve with length < 2.2.')

# Finally we might want to take a look at the tiling we've produced, so we should save
# it (as a .svg type image).
T.save_diagram('tutorial_tiling_diagram.svg')

# As the svg format only works with integer coordinates and we have a lot of fine detail
# we might actuall need to save it using higher resolution coordinates.
T.save_diagram('tutorial_tiling_diagram_high_res.svg', size=2000)


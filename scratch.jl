

using DataFrames
mam = readtable("/Users/michael/Google Drev/Mountain project/New Polygons/Species grid distributions/mammals_PA_matrix.csv")
mam[3] = map(x->"$x", mam[3])
mama = Assemblage(mam)

# these should be named vectors
ric = richness(mama)
occu = occupancy(mama)
nspecies(mama)
sitenames(mama) # This should be the actual sitenames - edit constructors?
maximum(richness(mama))

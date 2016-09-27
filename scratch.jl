
using DataFrames
mam = readtable("/Users/michael/Google Drev/Mountain project/New Polygons/Species grid distributions/mammals_PA_matrix.csv")
mam[3] = map(x->"$x", mam[3])


mama = Assemblage(mam)
## Problems with saving

import JLD
# Go to here to begin with
mamcop = JLD.load("mamobj.jld", "mamnod")

# these should be named vectors
ric = richness(mamcop)
occu = occupancy(mamcop)
nspecies(mamcop)
specnames(mamcop) # This should be the actual sitenames - edit constructors?
maximum(richness(mamcop))

println(createsummaryline(specnames(mamcop)))

mamcop

fieldnames(mamcop)

fieldnames(mamcop.site)

fieldnames(mamcop.occ)

#TODO need to write proper doc strings - use what is already in the R package

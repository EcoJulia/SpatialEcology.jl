

amp = readtable("data/amph_Europe.csv")
amp[3] = map(x->"$x", amp[3])

amp = Assemblage(amp)
## Problems with saving

plot(amp)
histogram(richness(amp))



using JLD

jldopen("mamobj.jld", "w") do file
   addrequire(file, SpatialEcology)
   write(file, "mam", mama)
end










# an easy way to get an object to work with
rdat = getRobject("/Users/michael/Documents/Projects/Current Projects/Core Corvoids/Data/corvids_Grid_nodiv.RData", "corvids")
corv = Assemblage(rdat)

using RCall
R"""
library(nodiv)
data(coquettes)
"""
coq = Assemblage(getRobject("coquettes"))

# 4-degree mammals
mam4 = getRobject("/Users/michael/Google Drev/Genetic diversity/Resubmission/terrestrial mammals data object.RData", "ter")
mam4 = Assemblage(mam4)


import JLD
# Go to here to begin with
mam = load("mamobj.jld", "mam")

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

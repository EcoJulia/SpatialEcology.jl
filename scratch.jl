
using DataFrames
mam = readtable("/Users/michael/Google Drive/Mountain project/New Polygons/Species grid distributions/mammals_PA_matrix.csv")
mam[3] = map(x->"$x", mam[3])


mama = Assemblage(mam)
## Problems with saving


using JLD

jldopen("mamobj.jld", "w") do file
   addrequire(file, SpatialEcology)
   write(file, "mam", mama)
end


# an easy way to get an object to work with
rdat = getRobject("/Users/michael/Documents/Projects/Current Projects/Core Corvoids/Data/corvids_Grid_nodiv.RData", "corvids")
corv = Assemblage(rdat)











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

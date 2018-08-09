using SpatialEcology
using DataFrames
using CSV
using Test


amphdat = CSV.read("../data/amph_Europe.csv")

# constructors
amph = Assemblage(amphdat[4:end], amphdat[1:3], sitecolumns = false)

# accesseors
@test extrema(richness(amph)) == (1, 20)
@test nsites(amph) == 1010
@test noccupied(amph) == 1010
@test nspecies(amph) == 73
@test noccurring(amph) == 73
@test occurring(amph, 718) == [15]
@test occupancy(amph)[1] == 353

# views
va = view(amph, species = 1:10)

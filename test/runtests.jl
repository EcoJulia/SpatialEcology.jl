using SpatialEcology
using DataFrames
using CSV
using Base.Test

# write your own tests here
amphdat = CSV.read("../data/amph_Europe.csv")
amph = Assemblage(amphdat[4:end], amphdat[1:3], sitecolumns = false)
@test extrema(richness(amph)) == (1, 20)
@test nsites(amph) == 1010
@test noccupied(amph) == 1010
@test nspecies(amph) == 73
@test noccurring(amph) == 73
@test occurring(amph, 718) == [15]

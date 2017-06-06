using SpatialEcology
using Base.Test

# write your own tests here
amphdat = readtable("../data/amph_Europe.csv")
amph = Assemblage(amphdat[4:end], amphdat[1:3])
@test extrema(richness(amph)) == (1, 20)
@test nsites(amph) == 1010
@test noccupied(amph) == 1010
@test nsites(amph) == 73
@test noccurring(amph) == 73
@test occurring(amph, 718) == [15]

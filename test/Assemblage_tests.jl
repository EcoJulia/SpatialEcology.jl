using DataFrames
using CSV
using SpatialEcology
using Test

@testset "Assemblage" begin
    amphdat = CSV.read(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv"), DataFrame)
    amph = Assemblage(amphdat[!,4:end], amphdat[!,1:3], sitecolumns = false)

    @test typeof(amph) == Assemblage{Bool,SpatialEcology.Locations{SpatialEcology.GridData}}

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

    #operations
    amp2 = aggregate(amph, 2)
    @test sum(richness(amp2)) == 2862
    @test nsites(amp2) == 285

    pointamph = Assemblage(amphdat[!,4:end], amphdat[!,1:3], sitecolumns = false, cdtype = SpatialEcology.pointdata)
    amp3 = aggregate(pointamph, 2)
    @test richness(amp3) == richness(amp2)
end

using SparseArrays
using Random
using SpatialEcology
using Distances
using StableRNGs
using Test

@testset "ComMatrix" begin
    rng = StableRNG(1337)

    datb = sprand(rng, Bool, 12,8,0.9)

    dati = sprand(rng, 8,13,0.9).*100
    dati = round.(Int, dati)

    datf = sprand(rng, 11,9,0.9)

    cmb = ComMatrix(datb)
    @test cmb isa ComMatrix{Bool}
    cmi = ComMatrix(dati,
        speciesnames = ["sp1", "sp2", "sp3", "sp4", "sp5", "sp6", "sp7", "sp8"],
        sitenames = [:a, :b, :c, :d, :e, :f, :g, :h, :i, :j, :k, :l, :m])
    @test cmi isa ComMatrix{<:Union{Int64, Int32}}
    cmf = ComMatrix(datf, sitecolumns = false)

    @test ComMatrix(Matrix(datf)).occurrences == ComMatrix(datf).occurrences

    @test occupancy(cmb) == [6, 6, 8, 7, 5, 8, 7, 8, 6, 7, 7, 8]
    @test occupancy(cmi) == [11, 13, 13, 10, 12, 12, 11, 11]
    @test occupancy(cmf) == [9, 10, 11, 11, 9, 9, 11, 9, 11]

    @test richness(cmb) == [9, 10, 12, 11, 11, 9, 11, 10]
    @test richness(cmi) == [8, 8, 7, 7, 7, 7, 6, 7, 8, 6, 8, 7, 7]
    @test richness(cmf) == [8, 9, 9, 8, 8, 7, 9, 8, 9, 8, 7]

    @test nsites(cmb) == 8
    @test nsites(cmi) == 13
    @test nsites(cmf) == 11

    @test nspecies(cmb) == 12
    @test nspecies(cmi) == 8
    @test nspecies(cmf) == 9

    @test occurring(cmb) == 1:12
    @test occurring(cmi) == 1:8
    @test occurring(cmf) == 1:9

    @test occupied(cmb) == 1:8
    @test occupied(cmi) == 1:13
    @test occupied(cmf) == 1:11

    @test noccurring(cmb) == 12
    @test noccurring(cmi) == 8
    @test noccurring(cmf) == 9

    @test noccupied(cmb) == 8
    @test noccupied(cmi) == 13
    @test noccupied(cmf) == 11

    @test getspecies(cmb, 3) == [true, true, true, true, true, true, true, true]
    @test getspecies(cmb, "species2") == [false, true, true, false, true, true, true, true]
    @test getspecies(cmi, 2) == [10, 61, 48, 6, 11, 69, 11, 36, 14, 91, 28, 93, 58]
    gcm = getspecies(cmf, 6)
    @test gcm[3] == 0.9426149574692801
    @test gcm isa SubArray

    @test speciesnames(cmi)[8] == "sp8"
    @test speciesnames(cmf)[3] == "species3"

    @test sitenames(cmi)[4] == "d"
    @test sitenames(cmb)[5] == "site5"

    @test length(sitenames(cmf)) == nsites(cmf)
    @test length(sitenames(cmi)) == nsites(cmi)
    @test length(sitenames(cmb)) == nsites(cmb)

    @test length(speciesnames(cmf)) == nspecies(cmf)
    @test length(speciesnames(cmi)) == nspecies(cmi)
    @test length(speciesnames(cmb)) == nspecies(cmb)

    @test sitetotals(cmb) == [9, 10, 12, 11, 11, 9, 11, 10]
    @test sitetotals(cmb) == richness(cmb)
    @test sitetotals(cmi) == [336, 520, 259, 301, 264, 243, 285, 450, 400, 280, 431, 470, 424]
    @test length(sitetotals(cmf)) == nsites(cmf)
    @test sitetotals(cmf)[5] ≈ 3.985356922989832

    @test speciestotals(cmb) == occupancy(cmb)
    @test speciestotals(cmi) == [572, 536, 703, 537, 729, 549, 514, 523]
    @test length(speciestotals(cmf)) == nspecies(cmf)
    @test speciestotals(cmf)[8] ≈ 5.171424049687359

    @test size(cmb) == (12, 8)
    @test size(cmi, 1) == 8

    @test cooccurring(cmb, 1, 3) == [true, false, true, true, true, false, true, true]
    @test cooccurring(cmf, [8, 3]) == [true, true, true, true, true, false, true, true, true, false, true]

    dist = pairwise(BrayCurtis(), view(cmi, sites=1:2))
    @test round.(dist, digits=1) == [0.0 0.3; 0.3 0.0]
  # getindex and setindex are to do

end

using Random
using Test
using SpatialEcology
using DataFrames
using SparseArrays

@testset "Views" begin
    Random.seed!(1337)

    datf = sprand(11,9,0.6)
    datf./=SpatialEcology.colsum(datf)
    com = ComMatrix(datf)

    vcom = view(com, species = [1,2,3,5,7,8], sites = [2,3,6,8,9])
    comc = copy(vcom)

    @test occupancy(vcom) == occupancy(comc)
    @test richness(vcom) == richness(comc)
    @test nsites(vcom) == nsites(comc)
    @test nspecies(vcom) == nspecies(comc)
    @test occurring(vcom) == occurring(comc)
    @test occupied(vcom) == occupied(comc)
    @test noccurring(vcom) == noccurring(comc)
    @test noccupied(vcom) == noccupied(comc)
    @test getspecies(vcom, 4) == getspecies(comc, 4)
    @test nrecords(vcom) == nrecords(comc)
    @test size(vcom) == size(comc)
    @test cooccurring(vcom, 1,3,4) == cooccurring(comc, [1,3,4])
    @test cooccurring(vcom, 3:5) == [false, true, false, true, false]

    @test SpatialEcology.asindices(1:5) == 1:5
    @test SpatialEcology.asindices([1,2,4]) == [1,2,4]
    @test SpatialEcology.asindices(1:3, 2:4) == 1:3
    @test SpatialEcology.asindices(["a", "c"], ["a", "b", "c"]) == [1,3]



end

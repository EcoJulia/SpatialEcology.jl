using SpatialEcology
using DataFrames
using CSV
using SparseArrays
using Test
import Random

Random.seed!(1337)
datb = sprand(Bool, 12,8,0.9)

dati = sprand(8,13,0.9).*100
dati = round.(Int, dati)

datf = sprand(11,9,0.9)


@testset "ComMatrix" begin
    cmb = ComMatrix(datb)
    @test cmb isa ComMatrix{Bool}
    cmi = ComMatrix(dati,
        specnames = ["sp1", "sp2", "sp3", "sp4", "sp5", "sp6", "sp7", "sp8"],
        sitenames = [:a, :b, :c, :d, :e, :f, :g, :h, :i, :j, :k, :l])
    @test cmi isa ComMatrix{<:Union{Int64, Int32}}
    cmf = ComMatrix(datf, sitecolumns = false)
    @test ComMatrix(Matrix(datf)).occurrences == ComMatrix(datf).occurrences

    @test occupancy(cmb) == [6, 6, 8, 7, 7, 8, 8, 8, 6, 6, 6, 7]
    @test occupancy(cmi) == [11, 13, 10, 12, 12, 13, 9, 13]
    @test occupancy(cmf) == [10, 10, 9, 10, 11, 8, 8, 10, 9]

    @test richness(cmb) == [12, 11, 10, 10, 11, 10, 12, 7]
    @test richness(cmi) == [7, 7, 8, 8, 5, 8, 8, 7, 8, 6, 6, 8, 7]
    @test richness(cmf) == [8, 8, 8, 9, 8, 8, 8, 7, 7, 7, 7]

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

    @test getspecies(cmb, 4) == [true, true, true, true, true, true, true, false]
    @test getspecies(cmi, 2) == [48, 89, 34, 39, 29, 73, 49, 23, 62, 71, 60, 57, 19]
    gcm = getspecies(cmf, 6)
    @test gcm[3] == 0.37980641711782925
    @test gcm isa SubArray

    @test specnames(cmi)[8] == "sp8"
    @test specnames(cmf)[3] == "species3"

    @test sitenames(cmi)[4] == "d"
    @test sitenames(cmb)[5] == "site5"

    @test length(sitenames(cmf)) == nsites(cmf)
    @test length(sitenames(cmi)) == nsites(cmi)
    @test length(sitenames(cmb)) == nsites(cmb)

    @test length(specnames(cmf)) == nspecies(cmf)
    @test length(specnames(cmi)) == nspecies(cmi)
    @test length(specnames(cmb)) == nspecies(cmb)






amphdat = CSV.read("../data/amph_Europe.csv")
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

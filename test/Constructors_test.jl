using DataFrames
using SparseArrays
using SpatialEcology
using Test

const SE = SpatialEcology

@testset "Constructors" begin

    @testset "ComMatrix from DataFrame" begin
        # Phylocom long format: [site, abundance, species] with a String species column
        phylo = DataFrame(site = ["s1", "s1", "s2", "s2"],
                          abu  = [1, 1, 1, 1],
                          species = ["a", "b", "a", "c"])
        cmp = ComMatrix(phylo)
        @test nspecies(cmp) == 3
        @test nsites(cmp) == 2
        @test sort(speciesnames(cmp)) == ["a", "b", "c"]
        @test sort(sitenames(cmp)) == ["s1", "s2"]

        # Wide format with species names in the first column
        wide = DataFrame(species = ["sp1", "sp2"], siteA = [1, 0], siteB = [1, 1])
        cm_s = ComMatrix(wide)
        @test speciesnames(cm_s) == ["sp1", "sp2"]
        @test sitenames(cm_s) == ["siteA", "siteB"]
        @test eltype(occurrences(cm_s)) == Bool

        # Wide format, no species column -> species get auto names, integer abundances
        cm_a = ComMatrix(DataFrame(siteA = [2, 0, 3], siteB = [0, 5, 1]))
        @test eltype(occurrences(cm_a)) <: Integer
        @test Matrix(occurrences(cm_a)) == [2 0; 0 5; 3 1]
        @test sitetotals(cm_a) == [5, 6]

        # Float abundances stay Float
        cm_r = ComMatrix(DataFrame(siteA = [0.2, 0.0, 0.3], siteB = [0.0, 0.5, 0.1]))
        @test eltype(occurrences(cm_r)) == Float64

        # Float values outside [0,1] still construct (with an informational note)
        cm_big = ComMatrix(DataFrame(siteA = [0.2, 1.5], siteB = [0.0, 2.0]))
        @test eltype(occurrences(cm_big)) == Float64
    end

    @testset "ComMatrix keyword constructor" begin
        cm = ComMatrix(species = ["a", "b", "a", "c"], sites = ["s1", "s1", "s2", "s2"])
        @test nspecies(cm) == 3
        @test nsites(cm) == 2

        cm2 = ComMatrix(species = ["a", "b", "a"], sites = ["s1", "s1", "s2"],
                        abundances = [3, 4, 5])
        @test eltype(occurrences(cm2)) <: Integer
        @test sum(occurrences(cm2)) == 12
    end

    @testset "ComMatrix from matrix with name validation" begin
        m = rand(Bool, 3, 4)
        cm = ComMatrix(m, ["a", "b", "c"], ["s1", "s2", "s3", "s4"])
        @test nspecies(cm) == 3 && nsites(cm) == 4
        @test_throws ArgumentError ComMatrix(m, ["a", "b"], ["s1", "s2", "s3", "s4"])
        @test_throws ArgumentError ComMatrix(m, ["a", "b", "c"], ["s1", "s2"])
    end

    @testset "Assemblage from a single worldmap DataFrame" begin
        wm = DataFrame(species = ["a", "a", "b", "c"],
                       ignored1 = [1, 2, 3, 4],
                       ignored2 = [5, 6, 7, 8],
                       long = [10.0, 11.0, 10.0, 12.0],
                       lat  = [55.0, 56.0, 55.0, 57.0])
        asm = Assemblage(wm)
        @test asm isa Assemblage
        @test nspecies(asm) == 3
        @test nsites(asm) == 3

        # a single DataFrame that is not in worldmap format is rejected
        @test_throws ErrorException Assemblage(DataFrame(a = 1:3, b = 4:6))
    end

    @testset "gridded coords with a gap" begin
        # x at 0,1,2,4 (the cell at x=3 is missing) exercises the grid-inference
        # branch that allows regularly spaced grids with empty cells
        coords = reduce(vcat, [Float64[x y] for x in [0.0, 1.0, 2.0, 4.0] for y in 0.0:1.0])
        cm = ComMatrix(rand(Bool, 3, size(coords, 1)), string.(1:3), string.(1:size(coords, 1)))
        a = Assemblage(cm, coords)
        @test a isa Assemblage{Bool, SE.Locations{SE.GridData}}
        @test xrange(a) == 0.0:1.0:4.0
    end

    @testset "dropemptyspecies via the (Locations, SpeciesData) constructor" begin
        gc = reduce(vcat, [Float64[x y] for x in 0.0:5.0 for y in 0.0:5.0])
        nsite = size(gc, 1)
        mat = zeros(Bool, 5, nsite); mat[1, 1] = true; mat[2, 2] = true   # species 3,4,5 empty
        cm = ComMatrix(mat, string.(1:5), string.(1:nsite))
        loc = SE.createLocations(gc, SE.griddata, DataFrame(sites = string.(1:nsite)))
        sp  = SE.SpeciesData(cm, DataFrame(name = string.(1:5)))
        a = Assemblage(loc, sp; dropemptyspecies = true)
        @test nspecies(a) == 2
    end

    @testset "Assemblage with coords as a DataFrame" begin
        cm = ComMatrix(rand(Bool, 4, 3), ["a", "b", "c", "d"], ["s1", "s2", "s3"])

        # two numeric columns -> taken directly as an [x y] matrix
        coords2 = DataFrame(x = [0.0, 1.0, 2.0], y = [0.0, 0.0, 0.0])
        a2 = Assemblage(cm, coords2)
        @test nsites(a2) == 3

        # three columns (a site-name column plus x and y) -> columns are guessed
        coords3 = DataFrame(site = ["s1", "s2", "s3"], x = [0.0, 1.0, 2.0], y = [0.0, 0.0, 0.0])
        a3 = Assemblage(cm, coords3)
        @test nsites(a3) == 3

        # an unusable shape is rejected
        bad = DataFrame(a = [1, 2, 3], b = [1, 2, 3], c = [1, 2, 3], d = [1, 2, 3])
        @test_throws ErrorException Assemblage(cm, bad)
    end

    @testset "testbool conversions" begin
        @test SE.testbool(1) === true
        @test SE.testbool(0) === false
        @test SE.testbool(true) === true
        @test SE.testbool(false) === false
        @test SE.testbool(missing) === false
        @test_throws ErrorException SE.testbool(2)       # integer > 1
        @test_throws ErrorException SE.testbool(1.5)     # non-integral number
        @test_throws ErrorException SE.testbool("x")     # non-numeric
    end

    @testset "inner-constructor dimension checks" begin
        cm = ComMatrix(rand(Bool, 3, 5), string.(1:3), string.(1:5))
        @test_throws DimensionMismatch SE.SpeciesData(cm, DataFrame(name = 1:2))
        coords = SE.GridData(reduce(vcat, [Float64[x y] for x in 0.0:1.0 for y in 0.0:1.0]))
        @test_throws DimensionMismatch SE.Locations{SE.GridData}(coords, DataFrame(sites = 1:99))
    end
end

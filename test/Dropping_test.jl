using SpatialEcology
using SparseArrays
using DataFrames
using Test

const SE = SpatialEcology
import SpatialEcology: places, things   # imported from EcoBase, not re-exported

@testset "Dropping empty species/sites and Locations indexing" begin
    # 6x6 grid (origin 0). Occurrences only for species 1,2,3 at sites 1,2,3,
    # so species 4,5,6 are empty and only 3 sites are occupied.
    gc = reduce(vcat, [Float64[x y] for x in 0.0:5.0 for y in 0.0:5.0])
    nsite = size(gc, 1)
    mat = zeros(Bool, 6, nsite)
    mat[1, 1] = true; mat[2, 2] = true; mat[3, 3] = true

    @testset "getindex(::Locations, inds)" begin
        for cdtype in (SE.griddata, SE.pointdata)
            asm = Assemblage(mat, gc, string.(1:nsite), string.(1:6); cdtype = cdtype)
            loc = places(asm)
            sub = loc[1:5]
            @test sub isa SE.Locations
            @test nsites(sub) == 5
            @test coordinates(sub) == coordinates(loc)[1:5, :]
        end
    end

    # the community matrix agrees with itself and with the traits and sitestats
    function aligned(asm)
        cm = commatrix(asm)
        speciesnames(cm) == traits(asm).name && sitenames(cm) == sitenames(asm) &&
            length(sitenames(cm)) == nsites(cm) &&
            SE.transposed(occurrences(cm)) == permutedims(SE._csc(occurrences(cm)))
    end

    @testset "dropempty* via the matrix constructor" begin
        a = Assemblage(mat, gc, string.(1:nsite), string.(1:6); dropemptyspecies = true)
        @test nspecies(a) == 3
        @test nsites(a) == nsite
        @test speciesnames(a) == speciesnames(commatrix(a)) == ["1", "2", "3"]
        @test aligned(a)

        b = Assemblage(mat, gc, string.(1:nsite), string.(1:6); dropemptysites = true)
        @test nsites(b) == 3
        @test nspecies(b) == 6
        @test size(coordinates(b)) == (3, 2)
        @test sitenames(commatrix(b)) == sitenames(b) == ["1", "2", "3"]
        @test aligned(b)

        c = Assemblage(mat, gc, string.(1:nsite), string.(1:6);
                       dropemptyspecies = true, dropemptysites = true)
        @test nspecies(c) == 3
        @test nsites(c) == 3
        @test aligned(c)
    end

    @testset "dropemptysites with an empty site between occupied ones" begin
        # s2 and s4 are empty, so s3 moves to position 2
        m = Bool[1 0 1 0; 0 0 1 0]
        pc = [0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0]
        for cdtype in (SE.griddata, SE.pointdata)
            asm = Assemblage(ComMatrix(m, ["a", "b"], ["s1", "s2", "s3", "s4"]), pc;
                             dropemptysites = true, cdtype = cdtype)
            cm = commatrix(asm)
            @test sitenames(cm) == sitenames(asm) == ["s1", "s3"]
            @test occupied(cm, "a") == [1, 2]
            @test occupied(cm, "b") == [2]
            @test richness(asm) == [1, 2]
            @test coordinates(asm) == pc[[1, 3], :]
            @test aligned(asm)
        end
    end

    @testset "dropemptysites works for grid and point data" begin
        for cdtype in (SE.griddata, SE.pointdata)
            a = Assemblage(mat, gc, string.(1:nsite), string.(1:6);
                           cdtype = cdtype, dropemptysites = true)
            @test nsites(a) == 3
            @test all(richness(a) .>= 1)            # every retained site is occupied
        end
    end

    @testset "dropemptysites via the (Locations, SpeciesData) constructor" begin
        loc = SE.createLocations(gc, SE.griddata, DataFrame(sites = string.(1:nsite)))
        sp  = SE.SpeciesData(ComMatrix(mat, string.(1:6), string.(1:nsite)),
                             DataFrame(name = string.(1:6)))
        a = Assemblage(loc, sp; dropemptysites = true)
        @test nsites(a) == 3
    end
end

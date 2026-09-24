using CSV
using DataFrames
using SpatialEcology
using SparseArrays
using Random
using StableRNGs
using Test

const SE = SpatialEcology
import SpatialEcology: places, things   # imported from EcoBase, not re-exported

@testset "Subsetting (views, copy, quantiles)" begin
    rng = StableRNG(55)
    gc = reduce(vcat, [Float64[x y] for x in 0.0:5.0 for y in 0.0:5.0])
    nsite = size(gc, 1)

    # an assemblage with deliberately empty species and empty sites
    mat = zeros(Bool, 6, nsite)
    mat[1, 1] = true; mat[2, 2] = true; mat[3, 3] = true   # species 4,5,6 and most sites empty
    asm = Assemblage(mat, gc, string.(1:nsite), string.(1:6))

    @testset "views with dropping" begin
        v = view(asm, dropempty = true)
        @test nspecies(v) == 3
        @test nsites(v) == 3

        @test nsites(view(asm, dropsites = true)) == 3
        @test nspecies(view(asm, dropspecies = true)) == 3

        # view of the SpeciesData component
        sv = view(things(asm))
        @test nspecies(sv) == nspecies(asm)
    end

    @testset "copy materialises every location type" begin
        for cdtype in (SE.griddata, SE.pointdata)
            a = Assemblage(Matrix(sprand(rng, Bool, 5, nsite, 0.5)), gc,
                           string.(1:nsite), string.(1:5); cdtype = cdtype)
            c = copy(a)
            @test typeof(c) == typeof(a)
            @test coordinates(c) == coordinates(a)
            @test occurrences(c) == occurrences(a)
        end
        @test copy(places(asm)) isa SE.Locations
        @test copy(commatrix(asm)) isa ComMatrix
        @test copy(things(asm)) isa SE.SpeciesData
    end

    @testset "Assemblage(view) materialises a SubAssemblage" begin
        v = view(asm, sites = 1:5, species = 1:4)
        a = Assemblage(v)
        @test a isa Assemblage
        @test nsites(a) == 5
        @test nspecies(a) == 4
    end

    @testset "asquantiles!" begin
        x = collect(1.0:20.0)
        SE.asquantiles!(x, 4)
        @test sort(unique(x)) == [1.0, 2.0, 3.0, 4.0]   # four quantile bins
    end

    # Names given as some AbstractString other than String (as CSV.jl reads
    # them) are stored as String, and selectors of any string type find them.
    # The inner constructor is used because the outer ones pass names through
    # `string`, which already turns SubStrings into Strings.
    @testset "views by name with non-String names" begin
        subnames(n) = [SubString("x$i", 2) for i in 1:n]
        a = Assemblage(ComMatrix{Bool}(sparse(mat), subnames(6), subnames(nsite)), gc)
        @test speciesnames(a) isa Vector{String}
        @test sitenames(a) isa Vector{String}

        @test speciesnames(view(a, species = subnames(2))) == ["1", "2"]
        v = view(a, species = ["2", "1"])
        @test speciesnames(v) == ["2", "1"]
        @test occurrences(v) == occurrences(asm)[[2, 1], :]
        v = view(a, sites = ["3", "1"])
        @test sitenames(v) == ["3", "1"]
        @test occurrences(v) == occurrences(asm)[:, [3, 1]]
        @test occupied(a, "2") == occupied(asm, "2")
        @test occurring(a, "3") == occurring(asm, "3")
    end
end

@testset "views by name on a CSV-read assemblage" begin
    occ = CSV.read(IOBuffer("site,abundance,species\ns1,1,sp_a\ns2,1,sp_a\ns2,1,sp_b\n"), DataFrame)
    coords = DataFrame(site = ["s1", "s2"], x = [0.0, 1.0], y = [0.0, 0.0])
    asm = Assemblage(occ, coords)
    @test speciesnames(asm) isa Vector{String}
    @test sitenames(asm) isa Vector{String}
    @test speciesnames(view(asm, species = ["sp_a"])) == ["sp_a"]
    # selecting with the CSV column itself, whose elements are not Strings
    @test speciesnames(view(asm, species = occ.species[3:3])) == ["sp_b"]
    @test sitenames(view(asm, sites = ["s2"])) == ["s2"]
    @test nspecies(view(asm, sites = ["s1"], dropspecies = true)) == 1
end

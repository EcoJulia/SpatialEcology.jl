using SpatialEcology
using SparseArrays
using DataFrames
using Random
using StableRNGs
using Test

const SE = SpatialEcology
import SpatialEcology: places, things   # imported from EcoBase, not re-exported

# Builds a fresh grid assemblage (mutating accessors run per-testset on a clean copy).
function makeasm(rng, nsp = 5)
    gc = reduce(vcat, [Float64[x y] for x in 0.0:5.0 for y in 0.0:5.0])
    nsite = size(gc, 1)
    mat = Matrix(sprand(rng, Bool, nsp, nsite, 0.5))
    Assemblage(mat, gc, string.(1:nsite), ["sp$i" for i in 1:nsp])
end

@testset "Get/Set data" begin
    rng = StableRNG(2024)

    @testset "accessors" begin
        asm = makeasm(rng)
        @test commatrix(asm) === commatrix(things(asm))
        @test occurrences(asm) === occurrences(commatrix(asm))
        @test traits(asm) === traits(things(asm))
        @test sitestats(asm) === places(asm).sitestats
        @test "name" in traitnames(asm)
        @test "sites" in sitestatnames(asm)
        @test size(coordinates(asm)) == (nsites(asm), 2)
    end

    @testset "dispersionfield" begin
        asm = makeasm(rng)
        df_idx = dispersionfield(asm, 1)
        @test length(df_idx) == nsites(asm)
        @test dispersionfield(asm, sitenames(asm)[1]) == df_idx   # index or name agree
    end

    @testset "addtraits! / addsitestats! with a vector" begin
        asm = makeasm(rng)
        addtraits!(asm, collect(1:nspecies(asm)), :score)
        @test traits(asm).score == collect(1:nspecies(asm))
        addsitestats!(asm, collect(1:nsites(asm)), :elev)
        @test sitestats(asm).elev == collect(1:nsites(asm))
        @test_throws ErrorException addtraits!(asm, [1, 2], :toofew)
        @test_throws ErrorException addsitestats!(asm, [1, 2], :toofew)
    end

    @testset "addtraits! / addsitestats! with a DataFrame join" begin
        asm = makeasm(rng)
        newtraits = DataFrame(name = speciesnames(asm), mass = collect(1.0:nspecies(asm)))
        addtraits!(asm, newtraits, :name)
        @test "mass" in traitnames(asm)
        @test collect(skipmissing(traits(asm).mass)) == collect(1.0:nspecies(asm))

        newsites = DataFrame(sites = sitenames(asm), elev = collect(1.0:nsites(asm)))
        addsitestats!(asm, newsites, :sites)
        @test "elev" in sitestatnames(asm)

        # no overlap in the join key -> abort
        @test_throws ErrorException addtraits!(asm, DataFrame(name = ["none"], mass = [0.0]), :name)
    end

    @testset "@traits and @sitestats" begin
        asm = makeasm(rng)
        addtraits!(asm, collect(1:nspecies(asm)), :score)
        addsitestats!(asm, collect(1:nsites(asm)), :elev)
        @test (@traits asm sum(:score)) == sum(1:nspecies(asm))
        @test (@sitestats asm sum(:elev)) == sum(1:nsites(asm))
        thr = 2                                            # a caller-local variable
        @test (@traits asm :score .> thr) == (collect(1:nspecies(asm)) .> thr)
    end

    @testset "getindex by symbol" begin
        asm = makeasm(rng)
        addsitestats!(asm, collect(1:nsites(asm)), :elev)
        addtraits!(asm, collect(1:nspecies(asm)), :score)
        @test asm[:elev] == collect(1:nsites(asm))         # a sitestat
        @test asm[:score] == collect(1:nspecies(asm))      # a trait
        @test asm[:, :elev] == asm[:elev]
        @test asm[!, :score] == asm[:score]      # `!` variant, trait branch
        @test asm[!, :elev] == asm[:elev]        # `!` variant, sitestat branch
        @test_throws ErrorException asm[:does_not_exist]
        @test_throws ErrorException asm[!, :does_not_exist]
    end
end

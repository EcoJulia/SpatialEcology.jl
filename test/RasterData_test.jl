using Test
using SpatialEcology
using Rasters
using DataFrames
using Random
using SpatialEcology: RasterData, SubRasterData

@testset "RasterData" begin
    rng = MersenneTwister(1)

    # --- a small 4x5 domain with some cells masked off ---
    x = X(1.0:4.0); y = Y(1.0:5.0)
    mask = Raster([(i + j) % 7 != 0 for i in 1:4, j in 1:5], (x, y))
    nsite = count(mask)

    # --- 3 species as Bool rasters over the same domain ---
    sp_names = ["sp_a", "sp_b", "sp_c"]
    ranges = map(1:3) do s
        Raster([(i * s + j) % 3 == 0 for i in 1:4, j in 1:5], (x, y); name = Symbol(sp_names[s]))
    end
    series = RasterSeries(ranges, (; name = sp_names))

    ss = DataFrame(pca1 = randn(rng, nsite), pca2 = randn(rng, nsite))
    asm = Assemblage(series, mask; sitestats = ss)

    @testset "construction" begin
        @test nspecies(asm) == 3
        @test nsites(asm) == nsite
        @test speciesnames(asm) == sp_names
        @test asm.site.coords isa RasterData
        @test size(coordinates(asm)) == (nsite, 2)
    end

    @testset "back-conversion" begin
        rr = richness_raster(asm)
        @test size(rr) == size(mask)
        # round-trip occurrences through a RasterSeries
        rs = to_rasterseries(asm)
        @test length(rs) == 3
        @test sum(sum(r[mask]) for r in rs) == sum(occurrences(asm))
    end

    @testset "site view" begin
        v = view(asm, sites = 1:5)
        @test v.site.coords isa SubRasterData
        @test nsites(v) == 5
        @test size(coordinates(v)) == (5, 2)
        @test size(occurrences(v)) == (3, 5)            # column-subview of parent
    end

    @testset "species + site view" begin
        v2 = view(asm, species = ["sp_a", "sp_c"], sites = 2:6)
        @test nspecies(v2) == 2
        @test nsites(v2) == 5
        @test speciesnames(v2) == ["sp_a", "sp_c"]
    end

    @testset "reordered view (a mask alone cannot represent this)" begin
        v3 = view(asm, sites = [5, 3, 1])
        @test coordinates(v3) == coordinates(asm)[[5, 3, 1], :]
    end

    @testset "subset to_raster keeps the full domain" begin
        v = view(asm, sites = 1:5)
        sub = to_raster(vec(sum(occurrences(v), dims = 1)), v)
        @test size(sub) == size(mask)
    end

    @testset "groupsites" begin
        addsitestats!(asm, repeat([:lo, :hi], outer = cld(nsite, 2))[1:nsite], :band)
        groups = groupsites(asm, :band)
        @test length(groups) == 2
        @test sum(nsites, groups) == nsite
    end

    @testset "copy(view) -> standalone RasterData" begin
        v3 = view(asm, sites = [5, 3, 1])
        c = copy(v3)
        @test c.site.coords isa RasterData
        @test coordinates(c) == coordinates(v3)
    end
end

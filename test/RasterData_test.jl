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

    @testset "dimension order (Y, X) is normalized" begin
        # Distinct X and Y ranges so a positional X/Y mix-up is detectable (and
        # would go out of bounds, since the grid is not square).
        xx = X(1.0:4.0); yy = Y(101.0:105.0)
        layers_xy = map(1:3) do s
            Raster([(i * s + j) % 3 == 0 for i in 1:4, j in 1:5], (xx, yy); name = Symbol(sp_names[s]))
        end
        mask_xy = Raster([(i + j) % 7 != 0 for i in 1:4, j in 1:5], (xx, yy))

        # the same data, transposed to (Y, X) order
        asm_yx = Assemblage(RasterSeries(permutedims.(layers_xy, Ref((Y, X))), (; name = sp_names)),
                            permutedims(mask_xy, (Y, X)))

        co = coordinates(asm_yx)
        @test all(in(1.0:4.0),     co[:, 1])    # X column holds X values, not Y
        @test all(in(101.0:105.0), co[:, 2])    # Y column holds Y values
        @test nsites(asm_yx) == nsite

        # identical set of sites to the (X, Y) build, regardless of input order
        asm_xy = Assemblage(RasterSeries(layers_xy, (; name = sp_names)), mask_xy)
        @test sortslices(co, dims = 1) == sortslices(coordinates(asm_xy), dims = 1)
    end

    @testset "RasterStack constructor" begin
        st = RasterStack((sp_a = ranges[1], sp_b = ranges[2], sp_c = ranges[3]))
        asm_st = Assemblage(st, mask)
        @test nspecies(asm_st) == 3
        @test speciesnames(asm_st) == sp_names
        @test occurrences(asm_st) == occurrences(asm)
    end

    @testset "default mask = full domain when no cells are missing" begin
        fullmask = Raster(trues(4, 5), (x, y))
        @test nsites(Assemblage(series)) == nsites(Assemblage(series, fullmask))
    end

    @testset "dims/extent mismatch is rejected" begin
        wrong = Raster([true for i in 1:3, j in 1:5], (X(1.0:3.0), Y(1.0:5.0)))
        @test_throws DimensionMismatch Assemblage(series, wrong)
    end

    @testset "aggregate" begin
        # qualified because Rasters also exports `aggregate` (name clash)
        agg = SpatialEcology.aggregate(asm, 2)
        @test agg isa Assemblage
        @test agg.site.coords isa RasterData
        @test nspecies(agg) == nspecies(asm)
        @test speciesnames(agg) == speciesnames(asm)
        @test nsites(agg) < nsites(asm)
        @test all(<=(nspecies(agg)), richness(agg))          # 'any' semantics
        # every occupied coarse cell must cover an occupied fine cell
        @test sum(richness(agg)) <= sum(richness(asm))
    end
end

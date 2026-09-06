using Test
using SpatialEcology
using Rasters
using DataFrames
using Random
using RecipesBase
using EcoBase
using SpatialEcology: RasterData, SubRasterData, indices, getcoords, places  # not exported

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

    @testset "domain / cellindices accessors" begin
        @test domain(asm) == Bool.(mask)
        @test domain(asm) === domain(asm.site) === domain(asm.site.coords)
        @test cellindices(asm) == findall(Bool.(mask))
        @test length(cellindices(asm)) == nsites(asm)
        # a value placed at site k lands in the k-th cell index
        vals = randn(rng, nsites(asm))
        r = to_raster(vals, asm)
        @test all(r[cellindices(asm)[k]] == vals[k] for k in 1:nsites(asm))
        # a view exposes its (reordered) sites but keeps the full parent domain
        v = view(asm, sites = [5, 3, 1])
        @test domain(v) == domain(asm)
        @test cellindices(v) == cellindices(asm)[[5, 3, 1]]
    end

    # indices() reads `cellinds`, not the `indices` field a GridData has, so a
    # RasterData needs its own method - without one every call here throws a
    # FieldError, including everything EcoBase derives from it.
    @testset "indices" begin
        @test size(indices(asm)) == (nsites(asm), 2)
        @test indices(asm, 1) == indices(asm)[:, 1]
        @test indices(asm, 2) == indices(asm)[:, 2]
        @test all(1 .<= indices(asm)[:, 1] .<= xcells(asm))
        @test all(1 .<= indices(asm)[:, 2] .<= ycells(asm))
        # the contract that ties them to coordinates(): index the ranges by
        # them and the coordinates come back, in the same site order
        @test xrange(asm)[indices(asm)[:, 1]] == coordinates(asm)[:, 1]
        @test yrange(asm)[indices(asm)[:, 2]] == coordinates(asm)[:, 2]
        # a reordered view carries its own site order through
        v = view(asm, sites = [5, 3, 1])
        @test indices(v) == indices(asm)[[5, 3, 1], :]
        @test xrange(v)[indices(v)[:, 1]] == coordinates(v)[:, 1]
    end

    # A north-up raster - Y ReverseOrdered, as anything read off disk with
    # Rasters normally is - must not be drawn upside down. The grid recipe
    # takes its y axis from yrange(grd, CellCentre()) but places cells by row
    # index, so the two have to run in the same direction; they only do
    # because xedges/yedges are read off the lookup in the lookup's own order
    # (see ext/RastersExt.jl). Derived edges ascend whatever the lookup does,
    # and pairing those with descending row indices flipped the map.
    @testset "north-up (descending Y) rasters are not drawn flipped" begin
        xd, yd = X(1.0:3.0), Y(3.0:-1.0:1.0)
        dom = Raster(trues(3, 3), (xd, yd))
        # one species, present only along y = 3.0 - the lookup's *first* row,
        # and the top of the map, so a flip is unmissable
        top = Raster([j == 1 for i in 1:3, j in 1:3], (xd, yd); name = :top)
        a = Assemblage(RasterSeries([top], (; name = ["top"])), dom)

        r = Float64.(richness(a))
        @test unique(coordinates(a)[:, 2][findall(==(1.0), r)]) == [3.0]

        # driven through apply_recipe so no plotting backend is needed, as in
        # PlotRecipes_test.jl
        _, ys, img = RecipesBase.apply_recipe(Dict{Symbol, Any}(), r,
                                              getcoords(places(a)))[1].args
        drawnrows = findall(row -> any(==(1.0), row),
                            collect(eachrow(replace(img, NaN => 0.0))))
        @test unique(collect(ys)[drawnrows]) == [3.0]
    end

    # The edges each grid reports, and the cell centres EcoBase derives from
    # them, follow the lookup rather than being rebuilt from a start and a
    # step - so they work for a descending axis and for a varying one.
    @testset "xedges/yedges follow the lookup" begin
        centres(e) = (e[1:(end - 1)] .+ e[2:end]) ./ 2
        for (label, xd, yd) in
            (("ascending", X(1.0:3.0), Y(1.0:3.0)),
             ("descending", X(1.0:3.0), Y(3.0:-1.0:1.0)),
             ("irregular", X([0.0, 1.0, 3.0]), Y([0.0, 2.0, 5.0])))
            rd = getcoords(places(Assemblage(
                RasterSeries([Raster([j == 1 for i in 1:3, j in 1:3], (xd, yd);
                                     name = :t)], (; name = ["t"])),
                Raster(trues(3, 3), (xd, yd)))))
            @test length(EcoBase.yedges(rd)) == ycells(rd) + 1
            @test length(EcoBase.xedges(rd)) == xcells(rd) + 1
            # the edges run the same way the grid's own range does
            @test issorted(EcoBase.yedges(rd)) == issorted(yrange(rd))
            # on a regular point lookup the centres are the coordinates back
            label == "irregular" ||
                @test centres(EcoBase.yedges(rd)) ≈ collect(yrange(rd))
        end
    end

    @testset "to_raster missingval keyword" begin
        counts = collect(1:nsites(asm))                  # an integer per-site vector
        default = to_raster(counts, asm)                 # historical behaviour
        @test eltype(default) == Int
        @test all(==(0), default[.!Bool.(mask)])         # off-domain cells are zero
        # NaN background promotes the output to float (off-domain renders missing)
        floated = to_raster(counts, asm; missingval = NaN)
        @test eltype(floated) <: AbstractFloat
        @test all(isnan, floated[.!Bool.(mask)])
        @test floated[cellindices(asm)] == Float64.(counts)
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

    @testset "coarsen" begin
        agg = coarsen(asm, 2)
        @test agg isa Assemblage
        @test agg.site.coords isa RasterData
        @test nspecies(agg) == nspecies(asm)
        @test speciesnames(agg) == speciesnames(asm)
        @test nsites(agg) < nsites(asm)
        @test all(<=(nspecies(agg)), richness(agg))          # 'any' semantics
        # every occupied coarse cell must cover an occupied fine cell
        @test sum(richness(agg)) <= sum(richness(asm))

        # Rasters' own `aggregate` verb is extended onto assemblages, so it works
        # unqualified here (SpatialEcology no longer exports `aggregate`) and
        # agrees with `coarsen`.
        @test occurrences(aggregate(asm, 2)) == occurrences(agg)
        @test occurrences(aggregate(maximum, asm, 2)) == occurrences(agg)
    end

    @testset "show" begin
        @test occursin("RasterData", sprint(show, asm.site.coords))
        @test occursin("SubRasterData view", sprint(show, view(asm, sites = 1:5).site.coords))
    end
end

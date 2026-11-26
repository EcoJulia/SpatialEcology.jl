using SpatialEcology
using SparseArrays
using Random
using StableRNGs
using Test

const SE = SpatialEcology
import SpatialEcology: places, things   # imported from EcoBase, not re-exported

# A 6x6 grid whose origin is at 0.0 coarsens cleanly by a factor of 2 into a 3x3 grid.
gridcoords() = reduce(vcat, [Float64[x y] for x in 0.0:5.0 for y in 0.0:5.0])

@testset "Operations (coarsen)" begin
    rng = StableRNG(99)
    gc = gridcoords()
    nsite = size(gc, 1)

    @testset "Boolean assemblage uses `any`" begin
        bmat = Matrix(sprand(rng, Bool, 6, nsite, 0.5))
        asm = Assemblage(bmat, gc, string.(1:nsite), string.(1:6))
        agg = coarsen(asm, 2)
        @test agg isa Assemblage{Bool}
        @test nsites(agg) == 9
        @test nspecies(agg) == 6
        # presence after coarsening is the OR of the grouped cells
        @test all(<=(1), occurrences(agg))
    end

    @testset "Integer assemblage: default (sum), explicit sum and a custom function" begin
        imat = round.(Int, rand(rng, 6, nsite) .* 4)
        asm = Assemblage(imat, gc, string.(1:nsite), string.(1:6))

        agg_default = coarsen(asm, 2)
        agg_sum     = coarsen(asm, 2, sum)
        agg_max     = coarsen(asm, 2, maximum)

        @test occurrences(agg_default) == occurrences(agg_sum)   # default for Int is `sum`
        @test sum(occurrences(agg_sum)) == sum(imat)             # abundance is conserved
        @test nsites(agg_sum) == 9
        @test all(occurrences(agg_max) .<= occurrences(agg_sum)) # max <= sum cell-by-cell
    end

    @testset "factor may be a tuple" begin
        bmat = Matrix(sprand(rng, Bool, 4, nsite, 0.5))
        asm = Assemblage(bmat, gc, string.(1:nsite), string.(1:4))
        @test nsites(coarsen(asm, (2, 2))) == nsites(coarsen(asm, 2))
    end

    @testset "coarsen Locations and GridTopology" begin
        bmat = Matrix(sprand(rng, Bool, 4, nsite, 0.5))
        asm = Assemblage(bmat, gc, string.(1:nsite), string.(1:4))

        lo = coarsen(places(asm), 2)
        @test lo isa SE.Locations{SE.GridData}
        @test nsites(lo) == 9

        gt = places(asm).coords.grid
        gt2 = coarsen(gt, 2)
        @test gt2 isa SE.GridTopology
        @test cells(gt2) == (3, 3)
    end

    points(occ, xy) = Assemblage(occ, xy, string.(axes(xy, 1)), ["a"]; cdtype = SE.pointdata)

    @testset "gridding point data" begin
        # the largest x and y fall on cell edges for most of these cellsizes
        pts = points([1 2 3 4], [0.5 0.5; 2.0 1.5; 1.0 0.25; 0.3 2.0])
        for cs in (1, 0.5, 0.1, (1.0, 0.5))
            gridded = coarsen(pts, cs)
            @test gridded isa Assemblage{Int, SE.Locations{SE.GridData}}
            @test SE.cellsize(places(gridded).coords.grid) == (first(cs), last(cs))
            @test sum(occurrences(gridded)) == 10
        end
        @test Matrix(occurrences(coarsen(pts, 1))) == [1 4 3 2]
        @test Matrix(occurrences(coarsen(pts, 2))) == [4 4 2]

        # all points share one x, on a cell edge
        @test nsites(coarsen(points([1 1], [1.0 0.0; 1.0 3.0]), 1)) == 2
        # the largest point is more than half a cell past the last full step from xmin
        @test sum(occurrences(coarsen(pts, 1; xmin = 0.3, ymin = 0.3))) == 10
        @test_throws ArgumentError coarsen(pts, 0)

        # rounded coordinates often sit on cell edges; each point gets a cell that holds it
        rounded = SE.PointData(round.(rand(rng, 50, 2) .* 20 .- 10, digits = 1))
        for cs in (0.1, 0.2, 0.5, 1)
            gt = SE._grid_from_factor(rounded, cs)
            ind = SE._apply_grid(rounded, gt)
            centres = hcat(gt.xs[ind[:, 1]], gt.ys[ind[:, 2]])
            @test all(abs.(rounded.coords .- centres) .<= cs / 2 + 1e-6)
        end
    end

    @testset "latitude 90" begin
        pts = points([1 1 1], [10.0 89.2; 10.5 90.0; 11.0 88.0])
        gridded = coarsen(pts, 1) # 89.2 and 90.0 share the top cell
        @test nsites(gridded) == 2
        @test last(places(gridded).coords.grid.ys) == 89.5
        # all points at the pole
        @test coordinates(coarsen(points([1 1], [0.0 90.0; 3.0 90.0]), 1))[:, 2] == [89.5, 89.5]
        # onto a given GridTopology
        gt = SE.GridTopology(range(10.5, step = 1.0, length = 2), range(88.5, step = 1.0, length = 2))
        @test nsites(coarsen(pts, gt)) == 2
        # these cannot be longitudes and latitudes, so 90 is an ordinary edge
        projected = coarsen(points([1 1], [500.0 89.2; 900.0 90.0]), 1)
        @test maximum(coordinates(projected)[:, 2]) == 90.5
    end

    @testset "gridded data needs an integer factor" begin
        asm = Assemblage(Matrix(trues(1, nsite)), gc, string.(1:nsite), ["a"])
        @test_throws ArgumentError coarsen(asm, 1.5)
        @test_throws ArgumentError coarsen(asm, 2.0)
        @test_throws ArgumentError coarsen(asm, (2, 1.5))
    end

    @testset "_default_fun" begin
        @test SE._default_fun(Bool) === any
        @test SE._default_fun(Int) === sum
        @test_throws ErrorException SE._default_fun(Float64)
    end
end

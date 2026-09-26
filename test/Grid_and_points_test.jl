using SpatialEcology
using Random
using SparseArrays
using StableRNGs
using Test

@testset "Grid" begin
    gr = SpatialEcology.creategrid([repeat(20.5 .+ (1:5), inner = 4) repeat(-30.5 .+ (2:2:8), outer = 5)])
    @test xmin(gr) == 21.5
    @test ymin(gr) == -28.5
    @test cellsize(gr) == (1.0, 2.0)
    @test cells(gr) == (5, 4)
    @test xrange(gr) == 21.5:1.0:25.5
    @test yrange(gr) == -28.5:2.0:-22.5
    @test boundingbox(gr).xmin == xmin(gr)
    @test sprint(show, boundingbox(gr)) == "xmin:\t21.5\nxmax:\t25.5\nymin:\t-28.5\nymax:\t-22.5\n"
end

@testset "Grids with cell sizes that are not exact in binary" begin
    function gridded(xs, ys; seed = 1)
        gc = reduce(vcat, [Float64[x y] for x in xs for y in ys])
        n = size(gc, 1)
        occ = Matrix(sprand(StableRNG(seed), Bool, 6, n, 0.1))
        Assemblage(occ, gc, string.(1:n), string.(1:6); cdtype = SpatialEcology.griddata), gc
    end

    grids = (("0.1 degree", range(-5.35, step = 0.1, length = 60), range(40.05, step = 0.1, length = 50)),
             ("1/12 degree", range(-20.0, step = 1/12, length = 60), range(30.0, step = 1/12, length = 50)),
             ("x 1.0, y 0.25", 1.0:1.0:30.0, range(0.125, step = 0.25, length = 40)))

    @testset "$desc" for (desc, xs, ys) in grids
        asm, gc = gridded(xs, ys)
        @test coordinates(asm) ≈ gc
        @test cells(asm) == (length(xs), length(ys))

        # copies of site views keep their coordinates
        rng = StableRNG(3)
        for _ in 1:200
            idx = sort(randperm(rng, nsites(asm))[1:rand(rng, 1:40)])
            c = copy(view(asm; sites = idx))
            @test coordinates(c) ≈ gc[idx, :]
        end

        # and so does dropping the empty sites in place
        dropped = copy(asm)
        keep = occupied(dropped)
        SpatialEcology.dropsites!(dropped.occ, dropped.site)
        @test coordinates(dropped) ≈ gc[keep, :]
    end
end

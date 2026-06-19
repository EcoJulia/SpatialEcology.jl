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

    @testset "_default_fun" begin
        @test SE._default_fun(Bool) === any
        @test SE._default_fun(Int) === sum
        @test_throws ErrorException SE._default_fun(Float64)
    end
end

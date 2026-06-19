using SpatialEcology
using SparseArrays
using Random
using StableRNGs
using Test

const SE = SpatialEcology
import SpatialEcology: places, things   # imported from EcoBase, not re-exported

@testset "show methods" begin
    rng = StableRNG(7)
    gc = reduce(vcat, [Float64[x y] for x in 0.0:5.0 for y in 0.0:5.0])
    nsite = size(gc, 1)
    asm = Assemblage(Matrix(sprand(rng, Bool, 6, nsite, 0.5)), gc, string.(1:nsite), string.(1:6))
    cm = commatrix(asm)

    @test occursin("ComMatrix with", sprint(show, cm))
    @test occursin("Assemblage with", sprint(show, asm))
    @test occursin("SubComMatrix with", sprint(show, view(cm, sites = 1:4)))
    @test occursin("SubAssemblage with", sprint(show, view(asm, sites = 1:4)))

    @test occursin("Spatial Grid", sprint(show, places(asm).coords))
    @test occursin("GridTopology", sprint(show, places(asm).coords.grid))
    @test occursin(r"xmin"i, sprint(show, boundingbox(asm)))

    addsitestats!(asm, repeat([:a, :b], outer = cld(nsite, 2))[1:nsite], :grp)
    groups = groupsites(asm, :grp)
    @test occursin("GroupedAssemblage with 2", sprint(show, groups))

    @testset "createsummaryline" begin
        @test SE.createsummaryline(["only"]) == "only"
        @test SE.createsummaryline(["a", "b", "c"]) == "a, b, c"
        @test occursin("...", SE.createsummaryline(["s$i" for i in 1:12]))
    end
end

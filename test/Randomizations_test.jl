using SpatialEcology
using SparseArrays
using Random
using StableRNGs
using Test

# Exercises src/Randomizations.jl, which previously had no test coverage at all.
# `curveball` randomization preserves both row sums (occupancy) and column sums
# (richness), so those marginals are the natural, seed-independent invariants.
@testset "Randomizations" begin
    rng = StableRNG(4242)

    @testset "Bool ComMatrix" begin
        cm = ComMatrix(sprand(rng, Bool, 15, 10, 0.4))
        gen = matrixrandomizer(cm)

        r = rand(gen)                       # returns a randomized copy
        @test r isa ComMatrix{Bool}
        @test sort(speciestotals(r)) == sort(speciestotals(cm))
        @test sort(sitetotals(r))    == sort(sitetotals(cm))

        r2 = rand!(gen)                     # mutates the generator's matrix in place
        @test r2 isa ComMatrix{Bool}
        @test sort(speciestotals(r2)) == sort(speciestotals(cm))
        @test sort(sitetotals(r2))    == sort(sitetotals(cm))
    end

    @testset "Bool Assemblage (grid and points)" begin
        gc = reduce(vcat, [Float64[x y] for x in 0.0:5.0 for y in 0.0:5.0])
        nsite = size(gc, 1)
        mat = Matrix(sprand(rng, Bool, 8, nsite, 0.4))

        for cdtype in (SpatialEcology.griddata, SpatialEcology.pointdata)
            asm = Assemblage(mat, gc, string.(1:nsite), string.(1:8); cdtype = cdtype)
            gen = matrixrandomizer(asm)

            a = rand(gen)
            @test a isa Assemblage{Bool}
            @test sort(richness(a))  == sort(richness(asm))
            @test sort(occupancy(a)) == sort(occupancy(asm))

            rand!(gen)                      # in-place variant should not throw
        end
    end

    @testset "non-Boolean falls back to a message" begin
        cmi = ComMatrix(round.(Int, sprand(rng, 6, 6, 0.6) .* 9 .+ 1))
        @test matrixrandomizer(cmi) isa AbstractString

        gc = reduce(vcat, [Float64[x y] for x in 0.0:3.0 for y in 0.0:3.0])
        ni = size(gc, 1)
        iasm = Assemblage(round.(Int, rand(rng, 4, ni) .* 4), gc, string.(1:ni), string.(1:4))
        @test matrixrandomizer(iasm) isa AbstractString
    end

    @testset "matrixrandomizations enum is available" begin
        @test matrixrandomizations isa DataType
        @test !isempty(instances(matrixrandomizations))
    end
end

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

    @testset "draws follow the generator's rng" begin
        cm = ComMatrix(sprand(StableRNG(1), Bool, 15, 10, 0.4))
        draw(seed) = rand(matrixrandomizer(cm, StableRNG(seed))).occurrences
        @test draw(7) == draw(7)
        @test draw(7) != draw(8)

        gc = reduce(vcat, [Float64[x y] for x in 0.0:5.0 for y in 0.0:5.0])
        nsite = size(gc, 1)
        mat = Matrix(sprand(StableRNG(2), Bool, 8, nsite, 0.4))
        asm = Assemblage(mat, gc, string.(1:nsite), string.(1:8))
        adraw(seed) = rand(matrixrandomizer(asm, StableRNG(seed))).occ.commatrix.occurrences
        @test adraw(7) == adraw(7)
        @test adraw(7) != adraw(8)

        # the global RNG is left alone
        Random.seed!(1)
        expected = rand()
        Random.seed!(1)
        draw(7)
        adraw(7)
        @test rand() == expected
    end

    @testset "species-wise queries see each draw" begin
        # the transposed occurrence matrix must be replaced along with the matrix
        cm = ComMatrix(sprand(StableRNG(3), Bool, 15, 10, 0.4))
        gen = matrixrandomizer(cm, StableRNG(5))
        for _ in 1:3
            r = rand!(gen)
            @test r.occurrences_t == permutedims(r.occurrences)
            @test SpatialEcology.thingoccurrences(r, 1) == r.occurrences[1, :]
        end
        @test rand!(gen) === rand!(gen)     # rand! reuses the generator's community
    end

    @testset "method = $method" for method in instances(matrixrandomizations)
        cm = ComMatrix(sprand(StableRNG(6), Bool, 8, 6, 0.4))
        r = rand(matrixrandomizer(cm, StableRNG(9); method))
        @test speciestotals(r) == speciestotals(cm)
        @test sitetotals(r)    == sitetotals(cm)
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

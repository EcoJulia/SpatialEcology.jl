
# Tests for view-of-view (nested subsetting) and related patterns that appear in the
# nodebased docs example. These were previously exercised only by the doctest, which
# meant breakage could go undetected between doc builds.
#
# Uses a small synthetic dataset to keep CI time negligible.

using SparseArrays, StableRNGs, Test, SpatialEcology
using Random: rand!

@testset "Nested views" begin
    rng = StableRNG(1337)

    # Small synthetic GridData-backed Assemblage: 5x4 grid (20 sites), 15 species
    xs = repeat(1.0:5.0, inner = 4)
    ys = repeat(1.0:4.0, outer = 5)
    coords = hcat(xs, ys)                          # 20×2, forms a regular grid
    occ = sprand(rng, Bool, 15, 20, 0.6)
    asm = Assemblage(ComMatrix(occ), coords)

    @test asm isa Assemblage{Bool, <:SpatialEcology.Locations{SpatialEcology.GridData}}

    sp_all  = speciesnames(asm)
    sp_top8 = sp_all[1:8]
    sp_top4 = sp_all[1:4]

    # --- First-level view ---
    sub1 = view(asm, species = sp_top8)
    @test sub1 isa SpatialEcology.SubAssemblage
    @test nspecies(sub1) == 8
    @test nsites(sub1)   == nsites(asm)

    # --- Nested view by species (the pattern that was broken) ---
    sub2 = view(sub1, species = sp_top4)
    @test sub2 isa SpatialEcology.SubAssemblage
    @test nspecies(sub2) == 4
    @test nsites(sub2)   == nsites(asm)
    @test occupancy(sub2) == occupancy(view(asm, species = sp_top4))
    @test richness(sub2)  == richness(view(asm, species = sp_top4))

    # --- Nested view by sites ---
    sub3 = view(sub1, sites = 1:10)
    @test nsites(sub3)   == 10
    @test nspecies(sub3) == 8

    # --- Nested view by both species and sites ---
    sub4 = view(sub1, species = sp_top4, sites = 1:10)
    @test nspecies(sub4) == 4
    @test nsites(sub4)   == 10

    # --- Three levels of nesting ---
    sub5 = view(sub4, sites = 1:5)
    @test nsites(sub5)   == 5
    @test nspecies(sub5) == 4

    # --- copy of a nested view produces a proper Assemblage with correct data ---
    cp = copy(sub2)
    @test cp isa Assemblage
    @test nspecies(cp)    == nspecies(sub2)
    @test nsites(cp)      == nsites(sub2)
    @test occupancy(cp)   == occupancy(sub2)
    @test richness(cp)    == richness(sub2)

    # --- matrixrandomizer on a SubAssemblage ---
    rmg = matrixrandomizer(sub1, StableRNG(99))
    randomized = rand(rmg)
    @test randomized isa Assemblage
    @test nspecies(randomized) == 8
    @test nsites(randomized)   == nsites(asm)
    # curveball preserves row and column sums
    @test occupancy(randomized) == occupancy(sub1)
    @test richness(randomized)  == richness(sub1)

    # --- view of a randomized Assemblage (the rand!(rmg) pattern from the docs) ---
    rand!(rmg)
    sub_of_rand = view(rmg.m, species = sp_top4)
    @test nspecies(sub_of_rand) == 4
end

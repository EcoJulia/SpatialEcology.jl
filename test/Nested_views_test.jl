
# Tests for view-of-view (nested subsetting) and related patterns that appear in the
# nodebased docs example. These were previously exercised only by the doctest, which
# meant breakage could go undetected between doc builds.

using CSV, DataFrames, SparseArrays, StableRNGs, Test, SpatialEcology
using Random: rand!

@testset "Nested views" begin
    amphdat = CSV.read(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv"), DataFrame)
    amph = Assemblage(amphdat[!,4:end], amphdat[!,1:3], sitecolumns = false)

    # --- First-level view ---
    sp_group1 = speciesnames(amph)[1:20]
    sp_group2 = speciesnames(amph)[1:10]

    sub1 = view(amph, species = sp_group1)
    @test sub1 isa SpatialEcology.SubAssemblage
    @test nspecies(sub1) == 20
    @test nsites(sub1) == nsites(amph)

    # --- Nested view: view of a SubAssemblage (the pattern that was broken) ---
    sub2 = view(sub1, species = sp_group2)
    @test sub2 isa SpatialEcology.SubAssemblage
    @test nspecies(sub2) == 10
    @test nsites(sub2) == nsites(amph)
    @test occupancy(sub2) == occupancy(view(amph, species = sp_group2))

    # --- View with site subsetting on a SubAssemblage ---
    sub3 = view(sub1, sites = 1:100)
    @test nsites(sub3) == 100
    @test nspecies(sub3) == 20

    # --- View with both species and sites on a SubAssemblage ---
    sub4 = view(sub1, species = sp_group2, sites = 1:50)
    @test nspecies(sub4) == 10
    @test nsites(sub4) == 50

    # --- copy of a nested view produces a proper Assemblage ---
    cp = copy(sub2)
    @test cp isa Assemblage
    @test nspecies(cp) == nspecies(sub2)
    @test nsites(cp) == nsites(sub2)
    @test occupancy(cp) == occupancy(sub2)

    # --- richness/occupancy of nested view matches the direct view ---
    @test richness(sub2) == richness(view(amph, species = sp_group2))
    @test occupancy(sub2) == occupancy(view(amph, species = sp_group2))

    # --- matrixrandomizer on a SubAssemblage ---
    rng = StableRNG(42)
    rand_sub = view(amph, species = sp_group1)
    rmg = matrixrandomizer(rand_sub, rng)
    randomized = rand(rmg)
    @test randomized isa Assemblage
    @test nspecies(randomized) == 20
    @test nsites(randomized) == nsites(amph)
    # curveball preserves row and column sums
    @test occupancy(randomized) == occupancy(rand_sub)
    @test richness(randomized) == richness(rand_sub)

    # --- view of a randomized (copy) Assemblage — the rand!(rmg) pattern ---
    rand!(rmg)
    sub_of_rand = view(rmg.m, species = sp_group2)
    @test nspecies(sub_of_rand) == 10

    # --- three levels of nesting (parent → clade → sister clade) ---
    parent = view(amph, species = sp_group1)
    child1 = view(parent, species = sp_group2)
    grandchild = view(child1, sites = 1:200)
    @test nspecies(grandchild) == 10
    @test nsites(grandchild) == 200
end

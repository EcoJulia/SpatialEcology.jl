using SpatialEcology
using SparseArrays
using Random
using StableRNGs
using RecipesBase
using Test

# The `plot(::Symbol, ::Assemblage)` recipe is exercised directly through
# RecipesBase.apply_recipe so the test does not need a plotting backend.
@testset "PlotRecipes" begin
    rng = StableRNG(11)
    gc = reduce(vcat, [Float64[x y] for x in 0.0:5.0 for y in 0.0:5.0])
    nsite = size(gc, 1)
    asm = Assemblage(Matrix(sprand(rng, Bool, 6, nsite, 0.5)), gc, string.(1:nsite), string.(1:6))
    addsitestats!(asm, Float64.(richness(asm)), :rich)

    res = RecipesBase.apply_recipe(Dict{Symbol, Any}(), :rich, asm)
    @test res isa Vector{RecipesBase.RecipeData}
    @test length(res) == 1
    @test res[1].args[1] == Float64.(richness(asm))

    # missing values in the plotted sitestat are turned into NaN
    addsitestats!(asm, [i == 1 ? missing : Float64(i) for i in 1:nsite], :withmiss)
    res2 = RecipesBase.apply_recipe(Dict{Symbol, Any}(), :withmiss, asm)
    vals = res2[1].args[1]
    @test isnan(vals[1])
    @test vals[2] == 2.0
end

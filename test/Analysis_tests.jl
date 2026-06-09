using DataFrames
using CSV
using SpatialEcology
using Test
using Plots

@testset "Analysis" begin
    amphdat = CSV.read(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv"), DataFrame)
    amph = Assemblage(amphdat[!,4:end], amphdat[!,1:3], sitecolumns = false)
    addtraits!(amph, asquantiles(occupancy(amph), 4), :quantile)
    gb = groupspecies(amph, :quantile) # error ?
    ps = [plot(g) for g in gb]
    p = ps[4];

    @test p[1][1][:seriestype] == :heatmap
    @test p[1][1][:z].surf[30,30] == 5.0

    # plot(vector, assemblage) recipe — the pattern used in the nodebased docs example
    # for plotting precomputed values (e.g. SOS) against site coordinates.
    vals = Float64.(richness(amph))
    pv = plot(vals, amph)
    @test pv[1][1][:seriestype] == :heatmap

    # also works on a SubAssemblage (e.g. rand_clade in the docs)
    sub = view(amph, species = 1:20)
    pv2 = plot(Float64.(richness(sub)), sub)
    @test pv2[1][1][:seriestype] == :heatmap
end

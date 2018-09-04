using DataFrames
using CSV
using SpatialEcology
using Test
using Plots

@testset "Analysis" begin
    amphdat = CSV.read(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv"))
    amph = Assemblage(amphdat[4:end], amphdat[1:3], sitecolumns = false)
    addtraits!(amph, asquantiles(occupancy(amph), 4), :quantile)
    gb = groupspecies(amph, :quantile)
    ps = [plot(g) for g in gb]
    p = ps[4];

    @test p[1][1][:seriestype] == :heatmap
    @test p[1][1][:z].surf[30,30] == 5.0
end

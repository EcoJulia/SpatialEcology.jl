using SpatialEcology
using Random
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

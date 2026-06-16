xmin(g::GridTopology) = minimum(g.xs)
ymin(g::GridTopology) = minimum(g.ys)
xcellsize(g::GridTopology) = step(g.xs)
ycellsize(g::GridTopology) = step(g.ys)
xcells(g::GridTopology) = length(g.xs)
ycells(g::GridTopology) = length(g.ys)
indices(g::SEGrid) = g.indices
indices(g::SEGrid, idx) = g.indices[:, idx]
boundingbox(g::GridTopology) = Bbox(xmin(g), xmax(g), ymin(g), ymax(g))

@forward_func GridData.grid xmin, ymin, xcellsize, ycellsize, xcells, ycells, boundingbox
@forward_func Locations{GridData}.coords xmin, ymin, xcellsize, ycellsize, xcells, ycells, boundingbox, indices
@forward_func SubGridData.grid xmin, ymin, xcellsize, ycellsize, xcells, ycells, boundingbox
@forward_func SubLocations{SubGridData}.coords xmin, ymin, xcellsize, ycellsize, xcells, ycells, boundingbox, indices
@forward_func SEAssemblage.site xmin, ymin, xcellsize, ycellsize, xcells, ycells, boundingbox, indices

# RasterData grid interface — derives topology from the mask raster's dims.
# Written once for both RasterData and SubRasterData (they share mask/cellinds).
# xrange/yrange return the lookup vectors; the others are scalars over the domain.
nplaces(rd::SubRasterData) = length(rd.cellinds)
EcoBase.xrange(rd::AnyRasterData) = parent(Rasters.lookup(rd.mask, Rasters.X()))
EcoBase.yrange(rd::AnyRasterData) = parent(Rasters.lookup(rd.mask, Rasters.Y()))
xmin(rd::AnyRasterData)       = minimum(EcoBase.xrange(rd))
xmax(rd::AnyRasterData)       = maximum(EcoBase.xrange(rd))
ymin(rd::AnyRasterData)       = minimum(EcoBase.yrange(rd))
ymax(rd::AnyRasterData)       = maximum(EcoBase.yrange(rd))
xcellsize(rd::AnyRasterData)  = abs(step(Rasters.lookup(rd.mask, Rasters.X())))
ycellsize(rd::AnyRasterData)  = abs(step(Rasters.lookup(rd.mask, Rasters.Y())))
xcells(rd::AnyRasterData)     = size(rd.mask, 1)
ycells(rd::AnyRasterData)     = size(rd.mask, 2)
boundingbox(rd::AnyRasterData) = Bbox(xmin(rd), xmax(rd), ymin(rd), ymax(rd))

@forward_func Locations{RasterData}.coords xmin, xmax, ymin, ymax, xcellsize, ycellsize, xcells, ycells, boundingbox
@forward_func SubLocations{SubRasterData}.coords xmin, xmax, ymin, ymax, xcellsize, ycellsize, xcells, ycells, boundingbox

show(io::IO, rd::RasterData) =
    println(io, "RasterData $(xcells(rd))×$(ycells(rd)) grid, $(nplaces(rd)) sites")
show(io::IO, rd::SubRasterData) =
    println(io, "SubRasterData view: $(nplaces(rd)) of $(count(rd.mask)) sites on a $(xcells(rd))×$(ycells(rd)) grid")


show(io::IO, b::Bbox) = println(io, "xmin:\t$(b.xmin)\nxmax:\t$(b.xmax)\nymin:\t$(b.ymin)\nymax:\t$(b.ymax)")
show(io::IO, g::GridData) = (println(io, "Spatial Grid\n", gridline(g.grid)); println(io, size(g.indices,1), " sites");)
show(io::IO, g::GridTopology) = println(io, "Spatial GridTopology\n", gridline(g))
gridline(g::GridTopology) =
    """
       lower left : $(xmin(g)), $(ymin(g))
       cellsizes  : $(xcellsize(g)), $(ycellsize(g))
       size       : $(xcells(g)), $(ycells(g))
    """

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


show(io::IO, b::Bbox) = println(io, "xmin:\t$(b.xmin)\nxmax:\t$(b.xmax)\nymin:\t$(b.ymin)\nymax:\t$(b.ymax)")
show(io::IO, g::GridData) = (println(io, "Spatial Grid\n", gridline(g.grid)); println(io, size(g.indices,1), " sites");)
show(io::IO, g::GridTopology) = println(io, "Spatial GridTopology\n", gridline(g))
gridline(g::GridTopology) =
    """
       lower left : $(xmin(g)), $(ymin(g))
       cellsizes  : $(xcellsize(g)), $(ycellsize(g))
       size       : $(xcells(g)), $(ycells(g))
    """

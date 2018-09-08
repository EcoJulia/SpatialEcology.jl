xmin(g::GridTopology) = g.xmin
ymin(g::GridTopology) = g.ymin
xcellsize(g::GridTopology) = g.xcellsize
ycellsize(g::GridTopology) = g.ycellsize
xcells(g::GridTopology) = g.xcells
ycells(g::GridTopology) = g.ycells
boundingbox(g::GridTopology) = Bbox(xmin(g), xmax(g), ymin(g), ymax(g))

@forward_func GridData.grid xmin, ymin, xcellsize, ycellsize, cellsize, xcells, ycells, cells, xrange, yrange, xmax, ymax, boundingbox
@forward_func Locations{GridData}.coords xmin, ymin, xcellsize, ycellsize, cellsize, xcells, ycells, cells, xrange, yrange, xmax, ymax, boundingbox


show(io::IO, b::Bbox) = println(io, "xmin:\t$(b.xmin)\nxmax:\t$(b.xmax)\nymin:\t$(b.ymin)\nymax:\t$(b.ymax)")
show(io::IO, g::GridData) = println(io,
    """
    Spatial grid
       lower left : $(xmin(g)), $(ymin(g))
       cellsizes  : $(xcellsize(g)), $(ycellsize(g))
       size       : $(xcells(g)), $(ycells(g))""")

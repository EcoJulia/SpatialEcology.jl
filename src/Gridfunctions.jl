xmin(g::GridTopology) = g.xmin
ymin(g::GridTopology) = g.ymin
cellsize(g::GridTopology) = g.xcellsize, g.ycellsize
xcellsize(g::GridTopology) = g.xcellsize
ycellsize(g::GridTopology) = g.ycellsize
xcells(g::GridTopology) = g.xcells
ycells(g::GridTopology) = g.ycells
cells(g::GridTopology) = g.xcells, g.ycells
xmax(g::GridTopology) = g.xmin + g.xcellsize*(g.xcells-1)
ymax(g::GridTopology) = g.ymin + g.ycellsize*(g.ycells-1)
xrange(g::GridTopology) = xmin(g):xcellsize(g):xmax(g)
yrange(g::GridTopology) = ymin(g):ycellsize(g):ymax(g)
boundingbox(g::GridTopology) = Bbox(xmin(g), xmax(g), ymin(g), ymax(g))
show(io::IO, b::Bbox) = println(io, "xmin:\t$(b.xmin)\nxmax:\t$(b.xmax)\nymin:\t$(b.ymin)\nymax:\t$(b.ymax)")

@forward_func EcoBase.AbstractGrid.grid xmin, ymin, cellsize, xcellsize, ycellsize, xcells, ycells, cells, xmax, ymax, xrange, yrange, boundingbox

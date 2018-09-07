xmin(g::GridTopology) = g.xmin
ymin(g::GridTopology) = g.ymin
xcellsize(g::GridTopology) = g.xcellsize
ycellsize(g::GridTopology) = g.ycellsize
xcells(g::GridTopology) = g.xcells
ycells(g::GridTopology) = g.ycells
boundingbox(g::GridTopology) = Bbox(xmin(g), xmax(g), ymin(g), ymax(g))

@forward_func SEGrid.grid xmin, ymin, xcellsize, ycellsize, cellsize, xcells, ycells, cells, xrange, yrange, xmax, ymax, boundingbox

show(io::IO, b::Bbox) = println(io, "xmin:\t$(b.xmin)\nxmax:\t$(b.xmax)\nymin:\t$(b.ymin)\nymax:\t$(b.ymax)")
show(io::IO, g::SEGrid) = println(io,
    """Spatial grid
       Lower left corner: $(xmin(g)), $(ymin(g))
       Cellsizes        : $(xcellsize(g)), $(ycellsize(g))
       Size             : $(xcells(g)), $(ycells(g))
       """)

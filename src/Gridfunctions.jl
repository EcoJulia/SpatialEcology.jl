xrange(g::GridTopology) = g.xmin:g.xcellsize:g.xmax
yrange(g::GridTopology) = g.ymin:g.ycellsize:g.ymax
bbox(g::GridTopology) = Bbox(g.xmin, g.xmin + g.xcells * g.xcellsize, g.ymin, g.xcells * g.xcellsize)
show(io::IO, b::Bbox) = println("xmin:\t$(g.xmin)\nxmax:\t$(g.xmax)\nymin:\t$(g.ymin)\nymax:\t$(g.ymax)\n")

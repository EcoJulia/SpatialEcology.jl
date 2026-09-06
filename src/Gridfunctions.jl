xmin(g::GridTopology) = minimum(g.xs)
ymin(g::GridTopology) = minimum(g.ys)
xcellsize(g::GridTopology) = step(g.xs)
ycellsize(g::GridTopology) = step(g.ys)
xcells(g::GridTopology) = length(g.xs)
ycells(g::GridTopology) = length(g.ys)
indices(g::SEGrid) = g.indices
indices(g::SEGrid, idx::Integer) = g.indices[:, idx]

# A RasterData holds its per-site cells as CartesianIndexes into the mask, not
# in an `indices` matrix, so the SEGrid methods above would reach for a field
# it does not have. The mask is always stored (X, Y), so a cell's
# CartesianIndex already reads as (x_index, y_index) and needs no reordering.
# This is the index counterpart of coordinates(rd::AnyRasterData), and takes
# the same route through `cellinds`, so the two stay in the same site order.
# It needs nothing from Rasters and so stays here rather than in the
# extension.
function indices(rd::AnyRasterData)
    return hcat([ci[1] for ci in rd.cellinds], [ci[2] for ci in rd.cellinds])
end
indices(rd::AnyRasterData, idx::Integer) = indices(rd)[:, idx]
boundingbox(g::GridTopology) = Bbox(xmin(g), xmax(g), ymin(g), ymax(g))

# xmin, ymin and indices are no longer forwarded. EcoBase 0.2 gives all three
# further methods - xmin(grd, anchor), indices(grd, order) - and @forward_func
# generates `f(x::T, args...)`, which is ambiguous with every one of them: the
# wrapper wins on argument one and EcoBase on argument two.
#
# Nothing replaces them above GridData, because nothing has to. EcoBase 0.2
# answers every gridded question for a type that holds gridded location data,
# so Locations, SubLocations and Assemblage get xmin, indices, xrange, xmax,
# the anchor forms and the ordered forms from the grid they hold, without
# this package restating any of it. Only the GridData -> GridTopology hop is
# invisible to EcoBase, so only that is written out.
@forward_func GridData.grid xcellsize, ycellsize, xcells, ycells, boundingbox
@forward_func Locations{GridData}.coords xcellsize, ycellsize, xcells, ycells, boundingbox
@forward_func SubGridData.grid xcellsize, ycellsize, xcells, ycells, boundingbox
@forward_func SubLocations{SubGridData}.coords xcellsize, ycellsize, xcells, ycells, boundingbox
@forward_func SEAssemblage.site xcellsize, ycellsize, xcells, ycells, boundingbox

xmin(g::GridData) = xmin(g.grid)
ymin(g::GridData) = ymin(g.grid)
xmin(g::SubGridData) = xmin(g.grid)
ymin(g::SubGridData) = ymin(g.grid)

# The RasterData grid interface (xrange/yrange/xmin/.../boundingbox, the
# Locations-level forwards, and the show methods) needs Rasters and lives in
# `ext/RastersExt.jl`. Only this Rasters-free method stays in core:
nplaces(rd::SubRasterData) = length(rd.cellinds)

show(io::IO, b::Bbox) = println(io, "xmin:\t$(b.xmin)\nxmax:\t$(b.xmax)\nymin:\t$(b.ymin)\nymax:\t$(b.ymax)")
show(io::IO, g::GridData) = (println(io, "Spatial Grid\n", gridline(g.grid)); println(io, size(g.indices,1), " sites");)
show(io::IO, g::GridTopology) = println(io, "Spatial GridTopology\n", gridline(g))
gridline(g::GridTopology) =
    """
       lower left : $(xmin(g)), $(ymin(g))
       cellsizes  : $(xcellsize(g)), $(ycellsize(g))
       size       : $(xcells(g)), $(ycells(g))
    """

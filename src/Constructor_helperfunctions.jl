
function isWorldmapData(dat::DataFrames.DataFrame, latlong = true)
  DataFrames.ncol(dat) == 5 || return false

  if eltypet(dat[:, 1]) <: AbstractString
    if eltypet(dat[:, 4]) <: Number
      if eltypet(dat[:, 5]) <: Number
        latlong || return true # risky
        if minimum(dropna(dat[:, 4])) > -181 && maximum(dropna(dat[:, 4])) < 181
          if minimum(dropna(dat[:, 5])) > -91 && maximum(dropna(dat[:, 5])) < 91
            return true
          end
        end
      end
    end
  end
  false
end


function parsesingleDataFrame(occ::DataFrames.DataFrame)
    if isWorldmapData(occ)
      println("Data format identified as Worldmap export file")
      coords = occ[4:5]
      coords[:sites] = createsitenames(coords)
      occ = DataFrame(site = coords[:sites], abu = ones(Int, DataFrames.nrow(occ)), species = occ[1])
      coords = unique(coords, :sites)
    else
        error("If coords is not specified, data must be in the worldmap format - of the column format [speciesnames, (ignored), (ignored), longitude, latitude]")
    end
  occ, coords
end

function guess_xycols(dat::DataFrames.DataFrame)
  numbers = map(x -> x<:Number, eltypest(dat))
  sum(.!numbers) == 1 || error("Site names cannot be numeric in the input matrix")
  ((findall(numbers)[1:2])...,)
end

function testbool(x::Int)
    x == 1 && return true
    x > 1 && error("Integer values > 1 used to create Boolean matrix")
    return false
end

testbool(x) = error("Value can not be interpreted as Boolean")
testbool(x::Missing) = false
testbool(x::Bool) = x
function testbool(x::Number)
  x == 0 && return false
  x == 1 && return true
  error("Value can not be interpreted as Boolean")
end

function dataFrametoSparseMatrix(dat::DataFrames.DataFrame, ::Type{T}) where T<:Bool
    is, js = Vector{Int}(), Vector{Int}()

    @inbounds for j in 1:DataFrames.ncol(dat)
        col = dat[:,j]
        for i in 1:DataFrames.nrow(dat)
            if testbool(col[i])
                push!(is, i)
                push!(js, j)
            end
        end
    end

    sparse(is, js, true, DataFrames.nrow(dat), DataFrames.ncol(dat))
end

function dataFrametoSparseMatrix(dat::DataFrames.DataFrame, ::Type{T}) where T<:Real
    is, js, vals = Vector{Int}(), Vector{Int}(), Vector{T}()

    @inbounds for j in 1:DataFrames.ncol(dat)
        col = dat[:,j]
        for i in 1:DataFrames.nrow(dat)
            if !ismissing(col[i]) && col[i] != 0
                push!(is, i)
                push!(js, j)
                push!(vals, col[i])
            end
        end
    end

    sparse(is, js, vals, DataFrames.nrow(dat), DataFrames.ncol(dat))
end

function match_commat_coords(occ::ComMatrix, coords::AbstractMatrix, sitestats::DataFrames.DataFrame)
  occ, coords, sitestats
 ## so far this does nothing TODO
end

function dropspecies!(occ::SpeciesData)
  occur = occurring(occ)
  occ.commatrix = occ.commatrix[occur, :]
  occ.traits = occ.traits[occur,:]
end

function dropbyindex!(site::Locations{PointData}, indicestokeep)
  site.coords = site.coords[indicestokeep,:]
  site.sitestats = site.sitestats[indicestokeep,:]
end

# these functions will be removed eventually
maxrange(x) = diff([extrema(x)...])[1]

# remember here - something wrong with the indices, make sure they are based from 1!

function dropbyindex!(site::Locations{GridData}, indicestokeep)
  site.indices = site.indices[indicestokeep,:]
  site.sitestats = site.sitestats[indicestokeep,:]
  site.grid.xmin = xrange(site.grid)[minimum(site.indices[:,1])]
  site.grid.ymin = yrange(site.grid)[minimum(site.indices[:,2])]
  site.grid.xcells = maxrange(site.indices[:,1]) + 1
  site.grid.ycells = maxrange(site.indices[:,2]) + 1
  site.indices = site.indices - minimum(site.indices) + 1
end

function dropsites!(occ::ComMatrix, site::SELocations)
  hasspecies = occupied(occ)
  occ.commatrix = occ.commatrix[:, hasspecies]
  dropbyindex!(site, hasspecies)
end

function createsitenames(coords::AbstractMatrix)
  size(coords, 2) == 2 || error("Only defined for matrices with two columns")
  mapslices(x->"$(x[1])_$(x[2])", coords, dims = 2)
end

function createsitenames(coords::DataFrames.DataFrame)
  size(coords, 2) == 2 || error("Only defined for matrices with two columns")
  ["$(coords[i,1])_$(coords[i,2])" for i in 1:DataFrames.nrow(coords)]
end

creategrid(coords::AbstractMatrix{<:Union{Number, Missing}}, tolerance = sqrt(eps())) =
    GridTopology(gridvar(coords[:,1], tolerance), gridvar(coords[:,2], tolerance))

# could allow for n-dimensional binning, using code from StatsBase.Histogram
function gridvar(x, tolerance = sqrt(eps()))
  sux = sort(unique(x))
  difx = diff(sux)
  length(difx) == 0 && error("Cannot make a grid with width 1 in the current implementation") #TODO
  rudifxmin, rudifxmax = extrema(unique(difx))
  if rudifxmax/rudifxmin < tolerance
    filter!(x->x>rudifxmax * tolerance, difx)
    rudifxmin, rudifxmax = extrema(unique(difx))
  end

  err1 = rudifxmax - rudifxmin
  if err1 > tolerance
    xx = [1, rudifxmax / rudifxmin]
    err2 = maximum(abs.(floor.(xx) - xx))
    err2  > tolerance && error("Cannot be converted to grid, as coordinate intervals are not constant. Try adjusting the tolerance (currently $tolerance)")
    difx = difx[difx .< rudifxmin + tolerance]
  end

  cellsize = mean(difx)
  min = minimum(sux)
  cellnumber = round(Int, maxrange(sux) / cellsize) + 1

  range(min, step = cellsize, length = cellnumber)
end

function getindices(coords::AbstractMatrix{<:Union{AbstractFloat, Missing}}, grid::GridTopology, tolerance = 2*sqrt(eps()))
  index1 = 1 .+ floor.(Int,(coords[:,1] .- xmin(grid)) ./ xcellsize(grid) .+ tolerance)
  index2 = 1 .+ floor.(Int,(coords[:,2] .- ymin(grid)) ./ ycellsize(grid) .+ tolerance)
  hcat(index1, index2)
end

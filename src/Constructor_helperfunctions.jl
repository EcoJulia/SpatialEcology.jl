

is01line{T <: Any}(vec::Union{AbstractVector{T},DataFrames.AbstractDataVector{T}}) = false
is01line{T <: Bool}(vec::Union{AbstractVector{T},DataFrames.AbstractDataVector{T}}) = true
is01line{T <: Number}(vec::Union{AbstractVector{T},DataFrames.AbstractDataVector{T}}) = length(setdiff(vec, [0, 1])) == 0


function BenHoltMatrix(commatrix::DataFrames.DataFrame)
  nc = DataFrames.ncol(commatrix)
  nc < 4 && return 0
  zeroonelines = vcat(DataFrames.colwise(is01line, commatrix)...)
  sum(zeroonelines) == 0 && return 0
  sum(zeroonelines[(nc-3):nc]) == 4 || return 0
  minimum(find(zeroonelines))
end


function isWorldmapData(dat::DataFrames.DataFrame, latlong = true)
  ncol(dat) == 5 || return false

  if eltype(dat[:, 1]) <: String
    if eltype(dat[:, 4]) <: Number
      if eltype(dat[:, 5]) <: Number
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
      if (firstnumeric = BenHoltMatrix(occ)) > 1
        println("Data assumed to be a concatenation of coordinates ($(firstnumeric - 1) columns) and occurrence matrix")
        coords = occ[1:(firstnumeric-1)]
        occ = occ[firstnumeric:end]
      else
        error("If not commatrix is already of type distrib_data or nodiv_data, a worldmap matrix, or a
          concatenation of coords and community matrix, coords must be specified")
      end
    end
  occ, coords
end

function parseDataFrame(occ::DataFrames.DataFrame)
  if ncol(occ) == 3 && eltypes(occ)[3] <: String
    println("Data format recognized as Phylocom")
    occ = unstack(occ, 1, 2)
  end

  if eltypes(occ)[1] <: String
    sites = Vector(occ[1])
    occ = occ[2:end]
  else
    sites = string.(1:DataFrames.nrow(occ))
  end

  try
    occ = dataFrametoNamedMatrix(occ, sites, Bool, dimnames = ("sites", "species"))
  catch
    occ = dataFrametoNamedMatrix(occ, sites, Int, dimnames = ("sites", "coordinates")) # This line means that this code is not completely type stable. So be it.
  end

  occ
end

function guess_xycols(dat::DataFrames.DataFrame)
  numbers = map(x -> x<:Number, eltypes(dat))
  sum(!numbers) == 1 || error("Site names cannot be numeric in the input matrix")
  ((find(numbers)[1:2])...)
end

function dataFrametoNamedMatrix(dat::DataFrames.DataFrame, rownames = string.(1:DataFrames.nrow(dat)), T::Type = Float64, replace = zero(T); sparsematrix = true, dimnames = ("A", "B"))
  colnames = string.(names(dat))
  a = 0
  for i in 1:ncol(dat)
    a += sum(isna(dat[i]))
    dat[i] = convert(Array, dat[i], replace)  #This takes out any NAs that may be in the data frame and replace with 0
  end

  a > 0 && println("$a NA values were replaced with $(replace)'s")
  try
    dat = sparsematrix ? sparse(Matrix{T}(dat)) : Matrix{T}(dat)
  catch
    error("Cannot convert DataFrame to Matrix{$T}")
  end

  dat = NamedArrays.NamedArray(dat, (Vector{String}(rownames), Vector{String}(colnames)), dimnames) #the vector conversion is a bit hacky
  dat
end



function match_commat_coords!(occ::ComMatrix, coords::AbstractMatrix{Float64}, sitestats::DataFrames.DataFrame)
  ()
 ## so far this does nothing TODO
end

function dropspecies!(occ::ComMatrix, traits::DataFrames.DataFrame)
  occurring = find(occupancy(occ) .> 0)
  occ = occ[:, occurring]
  traits = traits[occurring,:]
  occ, traits
end

function dropsites!(occ::ComMatrix, coords::AbstractMatrix, sitestats::DataFrames.DataFrame)
  hasspecies = find(richness(occ) .> 0)
  occ = occ[hasspecies,:]
  coords = coords[hasspecies,:]
  sitestats = sitestats[hasspecies,:]
  occ, coords, sitestats
end

function createsitenames(coords::AbstractMatrix)
  size(coords, 2) == 2 || error("Only defined for matrices with two columns")
  mapslices(x->"$(x[1])_$(x[2])", coords, 2)
end

function createsitenames(coords::DataFrames.DataFrame)
  size(coords, 2) == 2 || error("Only defined for matrices with two columns")
  ["$(coords[i,1])_$(coords[i,2])" for i in 1:DataFrames.nrow(coords)]
end


creategrid(coords::NamedArrays.NamedMatrix{Float64}, tolerance = sqrt(eps())) =
    GridTopology(gridvar(coords[:,1], tolerance)..., gridvar(coords[:,2], tolerance)...)


function gridvar(x, tolerance = sqrt(eps()))  # TODO this code is 'borrowed' from sp:::points2grid in R, which is GPL2, so cannot be distributed in this package
  sux = sort(unique(x))
  difx = diff(sux)
  length(difx) == 0 && error("Cannot make a grid with width 1 in the current implementation") #TODO
  rudifx = [extrema(unique(difx))...]
  if rudifx[1]/rudifx[2] < tolerance
    difx = difx[difx .> rudifx[2] * tolerance]
    rudifx = [extrema(unique(difx))...]
  end

  err1 = diff(rudifx)
  if err1 > tolerance
    xx = rudifx ./ minimum(rudifx)
    err2 = maximum(abs(floor(xx) - xx))
    err2  > tolerance && error("Cannot be converted to grid, as coordinate intervals are not constant. Try adjusting the tolerance (currently $tolerance)")
    difx = difx[difx .< rudifx[1] + tolerance]
  end

  cellsize = mean(difx)
  min = minimum(sux)
  cellnumber = Int(round(diff([extrema(sux)]) / cellsize) + 1)

  min, cellsize, cellnumber
end












isgrid(coords::AbstractMatrix) = isgridvar(coords[:,1]) && isgridvar(coords[:,2])

  function isgridvar(coord::AbstractVector)
    dists = diff(sort(unique(signif.(vec(coord), 8)))) # a bit hacky, prone to cause errors
    freqs = freq(dists)
    maxval = maximum(values(freqs))
    smallest = minimum(dists[dists .> 0])
    most_common = [k for (k,v) in freqs if v == maxval]
    sum(dists .% smallest) == 0 & length(intersect([smallest], most_common)) > 0
  end


  function freq{T}(v::AbstractVector{T})
    freqs = Dict{T, Int}()
    for i in v
      if haskey(freqs, i)
        freqs[i] += 1
      else
        freqs[i] = 0
      end
    end
    freqs
  end



is01line{T <: Any}(vec::AbstractDataVector{T}) = false
is01line{T <: Bool}(vec::AbstractDataVector{T}) = true
is01line{T <: Number}(vec::AbstractDataVector{T}) = length(setdiff(vec, [0, 1])) == 0


function BenHoltMatrix(commatrix::DataFrames.DataFrame)
  nc = ncol(commatrix)
  nc < 4 && return 0
  zeroonelines = vcat(colwise(is01line, commatrix)...)
  sum(zeroonelines) == 0 && return 0
  sum(zeroonelines[(nc-3):nc]) == 4 && return 0
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
      occ = DataFrame(site = coords[:sites], abu = ones(Int, nrow(occ)), species = occ[1])
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
  end
  occ, coords
end

function parseDataFrame(occ::DataFrames.DataFrame)
  if ncol(occ) == 3 && eltypes(occ)[3] <: String
    println("Data format recognized as Phylocom")
    tmp = unstack(occ, 1, 2)
  end

  if eltypes(occ)[1] <: String
    sites = Vector(occ[1])
    occ = occ[2:end]
  else
    sites = string.(1:nrow(occ))
  end

  for i in 1:ncol(occ)
    occ[i] = convert(Array, occ[i], 0)  #This takes out any NAs that may be in the data frame and replace with 0
  end

  try  # Let us see if it can be translated to a bool
      occ = Matrix{Bool}(tmp)
  catch
      occ = Matrix{Int}(tmp)  # This line means that this code is not completely type stable. So be it.
  end

  occ = NamedArrays.NamedArray(occ, sites, string.(names(occ)))
  occ
end



# constructor helper functions

function createNodeBySpeciesMatrix(tree::Phylogenetics.Phylogeny)
   colname = tree.tipLabel
   rowname = ["$num" for num in nodeNumbers(tree)]
   nodespecies = NamedArray(zeros(Int, nNode(tree), nTip(tree)), (rowname, colname))

   function loc(tree::Phylogenetics.Phylogeny, node::Int)
     if node <= nTip(tree)
       return node
     end
     ret = [loc(tree, descendants(node, tree)[1])..., loc(tree, descendants(node, tree)[2])...]
     nodespecies[node-nTip(tree), ret] = 1
     ret
   end
   loc(tree, basalNode(tree))

   nodespecies
end

function createsitenames(coords::AbstractMatrix)
  size(coords, 2) == 2 || error("Only defined for matrices with two columns")
  mapslices(x->"$(x[1])_$(x[2])", coords, 2)
end

function createsitenames(coords::DataFrames.DataFrame)
  size(coords, 2) == 2 || error("Only defined for matrices with two columns")
  ["$(coords[i,1])_$(coords[i,2])" for i in 1:nrow(coords)]
end

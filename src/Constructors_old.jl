

function OccMatrix{T <: Union{Int, Bool}}(x::AbstractMatrix{T})
    println("success")
    NamedArrays.NamedArray(sparse(x))
end

function OccMatrix{T <: Union{Int, Bool}}(x::AbstractMatrix{T}, species::Vector{String}, sites::Vector{String})
    NamedArrays.NamedArray(sparse(x), (sites, species))
end

# reconsidering this.
function OccMatrix{T <: Union{Int, Bool}}(x::Matrix{T}, species::Vector{String}, sites::Vector{String})
    NamedArrays.NamedArray(sparse(x), (sites, species))
end

#function PAMatrix(x::Matrix{Bool}; species = :auto, sites = :auto)
#    OccMatrix(x, species = species, sites = sites)
#end
#
#function AbundanceMatrix(x::Matrix{Int}; species = :auto, sites = :auto)
#    OccMatrix(x, species = species, sites = sites)
#end











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

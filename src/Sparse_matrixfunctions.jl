
# Generic column/row sums. These were historically provided by EcoBase, which no
# longer exports them; the Bool and SubArray specializations below override these
# fast paths. Shape (1×N / N×1) is kept to match the old EcoBase behaviour — call
# sites wrap these in `vec` (see sitetotals/speciestotals).
colsum(x::AbstractMatrix) = sum(x, dims = 1)
rowsum(x::AbstractMatrix) = sum(x, dims = 2)

# Functions for sparse array view sums - from https://discourse.julialang.org/t/slow-arithmetic-on-views-of-sparse-matrices/3644

function rowsum(x::SubArray{T,2,P}) where {T,P<:SparseMatrixCSC}
    pi = parentindices(x)
    (x.parent * sparse(pi[2], ones(Int, length(pi[2])), ones(Int, length(pi[2])), size(x.parent, 2), 1))[pi[1]]
end

function rowsum(x::SparseMatrixCSC{Bool})
    active = zeros(Int, size(x, 1))
    rv = rowvals(x)
    @inbounds for i in 1:nnz(x)
        active[rv[i]] += 1
    end
    active
end

colsum(x::SparseMatrixCSC{Bool}) = diff(SparseArrays.getcolptr(x))

function colsum(x::SubArray{T,2,P}) where {T,P<:SparseMatrixCSC}
    pi = parentindices(x)
    (sparse(ones(Int, length(pi[1])), pi[1], ones(Int, length(pi[1])), 1, size(x.parent, 1)) * x.parent)[pi[2]]
end

# Functions for finding nonzero rows and columns from Dan Getz, http://stackoverflow.com/questions/43968445/identify-which-rows-or-columns-have-values-in-sparse-matrix

function nzrows(a::SparseMatrixCSC)
    active = falses(size(a, 1))
    rv = rowvals(a)
    @inbounds for i in 1:nnz(a)
        active[rv[i]] = true
    end
    findall(active)
end

function nzcols(a::SparseMatrixCSC)
    cp = SparseArrays.getcolptr(a)
    findall(i -> cp[i] < cp[i+1], 1:size(a, 2))
end

inrange(v,r) = searchsortedlast(v,last(r))>=searchsortedfirst(v,first(r))

function sortedintersecting(v1, v2)
    i,j = firstindex(v1), firstindex(v2)
    @inbounds while i <= lastindex(v1) && j <= lastindex(v2)
        if v1[i] == v2[j] return true
        elseif v1[i] > v2[j] j += 1
        else i += 1
        end
    end
    return false
end

function nzcols(b::SubArray{T,2,P,Tuple{Vector{Int64},Vector{Int64}}} where {T,P<:SparseMatrixCSC}
  )
    pi = parentindices(b)
    brows = sort(unique(pi[1]))
    @inbounds return [k
      for (k,i) in enumerate(pi[2])
      if !isempty(nzrange(b.parent, i)) &&
        sortedintersecting(rowvals(b.parent)[nzrange(b.parent,i)], brows)]
end

function nzcols(b::SubArray{T,2,P,M} where {T,P<:SparseMatrixCSC, M<:Tuple{AbstractUnitRange{Int64},AbstractUnitRange{Int64}}}
  )
    pi = parentindices(b)
    @inbounds return collect(i+1-first(pi[2])
      for i in pi[2]
      if !isempty(nzrange(b.parent, i)) &&
        inrange(rowvals(b.parent)[nzrange(b.parent,i)], pi[1]))
end

function nzcols(b::SubArray{T,2,P,M} where {T,P<:SparseMatrixCSC, M<:Tuple{AbstractUnitRange{Int64},Vector{Int64}}}
  )
    pi = parentindices(b)
    @inbounds return [k
      for (k,i) in enumerate(pi[2])
      if !isempty(nzrange(b.parent, i)) &&
        inrange(rowvals(b.parent)[nzrange(b.parent,i)], pi[1])]
end

function nzcols(b::SubArray{T,2,P,M} where {T,P<:SparseMatrixCSC,M<:Tuple{Vector{Int64},AbstractUnitRange{Int64}}}
  )
    pi = parentindices(b)
    brows = sort(unique(pi[1]))
    @inbounds return collect(i+1-first(pi[2])
      for i in pi[2]
      if !isempty(nzrange(b.parent, i)) &&
        sortedintersecting(rowvals(b.parent)[nzrange(b.parent,i)], brows))
end

function findin2(inds,v,w)
    i, j = firstindex(v), firstindex(w)
    res = Vector{Int}()
    @inbounds while i<=lastindex(v) && j<=lastindex(w)
        if v[i]==w[j]
            push!(res,inds[i])
            i += 1
        elseif (v[i]<w[j]) i += 1
        else j += 1
        end
    end
    return res
end

function nzrows(b::SubArray{T,2,P,Tuple{Vector{Int64}, U}} where {T,P<:SparseMatrixCSC, U<:Union{M, Vector{Int}}} where {M <: AbstractUnitRange{Int64}}
  )
    pi = parentindices(b)
    active = falses(length(pi[1]))
    inds = sortperm(pi[1])
    brows = (pi[1])[inds]
    @inbounds for c in pi[2]
      active[findin2(inds, brows, rowvals(b.parent)[nzrange(b.parent, c)])] .= true
    end
    return findall(active)
end

function nzrows(b::SubArray{T,2,P,Tuple{M, U}} where {T,P<:SparseMatrixCSC, U<:Union{M, Vector{Int}}} where M <: AbstractUnitRange{Int64}
  )
    pi = parentindices(b)
    rv = rowvals(b.parent)
    active = falses(length(pi[1]))
    @inbounds for c in pi[2]
        for r in nzrange(b.parent, c)
            if rv[r] in pi[1]
                active[rv[r]+1-first(pi[1])] = true
            end
        end
    end
    return findall(active)
end

_nnz(x::AbstractMatrix) = nnz(x)
function _nnz(b::SubArray{T,2,P}) where {T,P<:SparseMatrixCSC}
    pi = parentindices(b)
    rv = rowvals(b.parent)
    ret = 0
    @inbounds for c in pi[2]
        for r in nzrange(b.parent, c)
            if rv[r] in pi[1]
                ret += 1
            end
        end
    end
    return ret
end

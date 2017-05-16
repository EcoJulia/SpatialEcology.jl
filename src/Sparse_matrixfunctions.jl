
# Functions for sparse array view sums - from https://discourse.julialang.org/t/slow-arithmetic-on-views-of-sparse-matrices/3644

rowsum(x) = sum(x,2)
rowsum{T,P<:SparseMatrixCSC}(x::SubArray{T,2,P}) = (x.parent * sparse(x.indexes[2],ones(Int,length(x.indexes[2])), ones(Int,length(x.indexes[2])), size(x.parent,2),1))[x.indexes[1]]
colsum(x) = sum(x,1)
colsum{T,P<:SparseMatrixCSC}(x::SubArray{T,2,P}) = (sparse(ones(Int,length(x.indexes[1])), x.indexes[1], ones(Int,length(x.indexes[1])),1,size(x.parent,1))*x.parent)[x.indexes[2]]

# Functions for finding nonzero rows and columns from Dan Getz, http://stackoverflow.com/questions/43968445/identify-which-rows-or-columns-have-values-in-sparse-matrix

function nzrows(a::SparseMatrixCSC)
    active = falses(a.m)
    for r in a.rowval
        active[r] = true
    end
    return find(active)
end

function nzcols(a::SparseMatrixCSC)
    res=Vector{Int}()
    foldl((x,y)->(if (x<a.colptr[y]) push!(res,y-1) end; a.colptr[y]),a.colptr[1],2:a.n+1)
    res
end

inrange(v,r) = searchsortedlast(v,last(r))>=searchsortedfirst(v,first(r))

function sortedintersecting(v1, v2)
    i,j = start(v1), start(v2)
    while i <= length(v1) && j <= length(v2)
        if v1[i] == v2[j] return true
        elseif v1[i] > v2[j] j += 1
        else i += 1
        end
    end
    return false
end

function nzcols{T,P<:SparseMatrixCSC}(b::SubArray{T,2,P,Tuple{Vector{Int64},Vector{Int64}}}
  )
    brows = sort(unique(b.indexes[1]))
    return [k
      for (k,i) in enumerate(b.indexes[2])
      if b.parent.colptr[i]<b.parent.colptr[i+1] &&
        sortedintersecting(b.parent.rowval[nzrange(b.parent,i)], brows)]
end

function nzcols{T,P<:SparseMatrixCSC}(b::SubArray{T,2,P,Tuple{UnitRange{Int64},UnitRange{Int64}}}
  )
    return collect(i+1-start(b.indexes[2])
      for i in b.indexes[2]
      if b.parent.colptr[i]<b.parent.colptr[i+1] &&
        inrange(b.parent.rowval[nzrange(b.parent,i)], b.indexes[1]))
end

function nzcols{T,P<:SparseMatrixCSC}(b::SubArray{T,2,P,Tuple{UnitRange{Int64},Vector{Int64}}}
  )
  return [k
    for (k,i) in enumerate(b.indexes[2])
    if b.parent.colptr[i]<b.parent.colptr[i+1] &&
        inrange(b.parent.rowval[nzrange(b.parent,i)], b.indexes[1])]
end

function nzcols{T,P<:SparseMatrixCSC}(b::SubArray{T,2,P,Tuple{Vector{Int64},UnitRange{Int64}}}
  )
    brows = sort(unique(b.indexes[1]))
    return collect(i+1-start(b.indexes[2])
      for i in b.indexes[2]
      if b.parent.colptr[i]<b.parent.colptr[i+1] &&
        sortedintersecting(b.parent.rowval[nzrange(b.parent,i)], brows))
end

function findin2(inds,v,w)
    i,j = start(v),start(w)
    res = Vector{Int}()
    while i<=length(v) && j<=length(w)
        if v[i]==w[j]
            push!(res,inds[i])
            i += 1
        elseif (v[i]<w[j]) i += 1
        else j += 1
        end
    end
    return res
end

function nzrows{T,P<:SparseMatrixCSC, U<:Union{UnitRange{Int64}, Vector{Int}}}(b::SubArray{T,2,P,Tuple{Vector{Int64}, U}}
  )
    active = falses(length(b.indexes[1]))
    inds = sortperm(b.indexes[1])
    brows = (b.indexes[1])[inds]
    for c in b.indexes[2]
      active[findin2(inds,brows,b.parent.rowval[nzrange(b.parent,c)])] = true
    end
    return find(active)
end

function nzrows{T,P<:SparseMatrixCSC, U<:Union{UnitRange{Int64}, Vector{Int}}}(b::SubArray{T,2,P,Tuple{UnitRange{Int64}, U}}
  )
    active = falses(length(b.indexes[1]))
    for c in b.indexes[2]
        for r in nzrange(b.parent,c)
            if b.parent.rowval[r] in b.indexes[1]
                active[b.parent.rowval[r]+1-start(b.indexes[1])] = true
            end
        end
    end
    return find(active)
end

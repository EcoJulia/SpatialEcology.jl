
# Functions for sparse array view sums - from https://discourse.julialang.org/t/slow-arithmetic-on-views-of-sparse-matrices/3644

rowsum(x) = sum(x,dims=2)
rowsum(x::SubArray{T,2,P}) where {T,P<:SparseMatrixCSC} = (x.parent * sparse(x.indexes[2],ones(Int,length(x.indexes[2])), ones(Int,length(x.indexes[2])), size(x.parent,2),1))[x.indexes[1]]
colsum(x) = sum(x,dims=1)
colsum(x::SubArray{T,2,P}) where {T,P<:SparseMatrixCSC} = (sparse(ones(Int,length(x.indexes[1])), x.indexes[1], ones(Int,length(x.indexes[1])),1,size(x.parent,1))*x.parent)[x.indexes[2]]

# Functions for finding nonzero rows and columns from Dan Getz, http://stackoverflow.com/questions/43968445/identify-which-rows-or-columns-have-values-in-sparse-matrix

function nzrows(a::SparseMatrixCSC)
    active = falses(a.m)
    @inbounds for r in a.rowval
        active[r] = true
    end
    return findall(active)
end

function nzcols(a::SparseMatrixCSC)
    res=Vector{Int}()
    foldl((x,y)->(if (x<a.colptr[y]) push!(res,y-1) end; a.colptr[y]),2:a.n+1, init=a.colptr[1])
    res
end

inrange(v,r) = searchsortedlast(v,last(r))>=searchsortedfirst(v,first(r))

function sortedintersecting(v1, v2)
    i,j = start(v1), start(v2)
    @inbounds while i <= length(v1) && j <= length(v2)
        if v1[i] == v2[j] return true
        elseif v1[i] > v2[j] j += 1
        else i += 1
        end
    end
    return false
end

function nzcols(b::SubArray{T,2,P,Tuple{Vector{Int64},Vector{Int64}}} where {T,P<:SparseMatrixCSC}
  )
    brows = sort(unique(b.indexes[1]))
    @inbounds return [k
      for (k,i) in enumerate(b.indexes[2])
      if b.parent.colptr[i]<b.parent.colptr[i+1] &&
        sortedintersecting(b.parent.rowval[nzrange(b.parent,i)], brows)]
end

function nzcols(b::SubArray{T,2,P,Tuple{UnitRange{Int64},UnitRange{Int64}}} where {T,P<:SparseMatrixCSC}
  )
    @inbounds return collect(i+1-start(b.indexes[2])
      for i in b.indexes[2]
      if b.parent.colptr[i]<b.parent.colptr[i+1] &&
        inrange(b.parent.rowval[nzrange(b.parent,i)], b.indexes[1]))
end

function nzcols(b::SubArray{T,2,P,Tuple{UnitRange{Int64},Vector{Int64}}} where {T,P<:SparseMatrixCSC}
  )
  @inbounds return [k
    for (k,i) in enumerate(b.indexes[2])
    if b.parent.colptr[i]<b.parent.colptr[i+1] &&
        inrange(b.parent.rowval[nzrange(b.parent,i)], b.indexes[1])]
end

function nzcols(b::SubArray{T,2,P,Tuple{Vector{Int64},UnitRange{Int64}}} where {T,P<:SparseMatrixCSC}
  )
    brows = sort(unique(b.indexes[1]))
    @inbounds return collect(i+1-start(b.indexes[2])
      for i in b.indexes[2]
      if b.parent.colptr[i]<b.parent.colptr[i+1] &&
        sortedintersecting(b.parent.rowval[nzrange(b.parent,i)], brows))
end

function findin2(inds,v,w)
    i,j = start(v),start(w)
    res = Vector{Int}()
    @inbounds while i<=length(v) && j<=length(w)
        if v[i]==w[j]
            push!(res,inds[i])
            i += 1
        elseif (v[i]<w[j]) i += 1
        else j += 1
        end
    end
    return res
end

function nzrows(b::SubArray{T,2,P,Tuple{Vector{Int64}, U}} where {T,P<:SparseMatrixCSC, U<:Union{UnitRange{Int64}, Vector{Int}}}
  )
    active = falses(length(b.indexes[1]))
    inds = sortperm(b.indexes[1])
    brows = (b.indexes[1])[inds]
    @inbounds for c in b.indexes[2]
      active[findin2(inds,brows,b.parent.rowval[nzrange(b.parent,c)])] = true
    end
    return findall(active)
end

function nzrows(b::SubArray{T,2,P,Tuple{UnitRange{Int64}, U}} where {T,P<:SparseMatrixCSC, U<:Union{UnitRange{Int64}, Vector{Int}}}
  )
    active = falses(length(b.indexes[1]))
    @inbounds for c in b.indexes[2]
        for r in nzrange(b.parent,c)
            if b.parent.rowval[r] in b.indexes[1]
                active[b.parent.rowval[r]+1-start(b.indexes[1])] = true
            end
        end
    end
    return findall(active)
end

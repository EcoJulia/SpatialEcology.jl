"""
    asquantiles!(x, n)
"""
function asquantiles!(x::AbstractVector, n::Int)
    quants = nquantile(x, n)[1:end-1]
    for (i,j) in enumerate(x)
        x[i] = searchsortedlast(quants, j)
    end
end

"""
    asquantiles(x, n)
"""
function asquantiles(x::AbstractVector, n::Int)
    quants = nquantile(x, n)[1:end-1]
    ret = similar(x)
    for (i,j) in enumerate(x)
        ret[i] = searchsortedlast(quants, j)
    end
    ret
end

eltypest(df) = map(x->Base.nonmissingtype(x), eltype.(eachcol(df)))
eltypet(x) = Base.nonmissingtype(eltype(x))

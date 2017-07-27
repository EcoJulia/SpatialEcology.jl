
function asquantiles!(x::AbstractVector, n::Int)
    quants = nquantile(x, n)
    for (i,j) in enumerate(x)
        x[i] = searchsortedlast(quants, j)
    end
end

function asquantiles(x::AbstractVector, n::Int)
    quants = nquantile(x, n)
    ret = similar(x)
    for (i,j) in enumerate(x)
        ret[i] = searchsortedlast(quants, j)
    end
    ret
end

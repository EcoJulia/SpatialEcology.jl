Tables.istable(::Type{<:AbstractComMatrix}) = true
Tables.rowaccess(::Type{<:AbstractComMatrix}) = true
Tables.columnaccess(::Type{<:AbstractComMatrix}) = true

struct ComMatrixRow{T}
    x::T
    names::Vector{Symbol}
end

# # Support indexing by rows
# struct AbstractComMatrixTable{transposed, T}
#     matrix::T
# end

AbstractComMatrixTable(acm::AbstractComMatrix; transpose::Bool=false) = AbstractComMatrixTable{transpose, typeof(acm)}(acm)

Base.propertynames(c::ComMatrixRow) = getfield(c, :names)
Base.getproperty(c::ComMatrixRow, nm::Symbol) = getfield(c, :x)[findfirst(==(nm), getfield(c, :names))]

function Base.iterate(acm::AbstractComMatrix, st=1)
    st > size(acm, 1) || return nothing
    return ComMatrixRow(acm[st, :], sitenames(acm)), st + 1
end

Base.IteratorEltype(::Type{AbstractComMatrix}) = Base.HasEltype()
Base.eltype(acm::AbstractComMatrix) = ComMatrixRow
Base.IteratorSize(::Type{AbstractComMatrix}) = Base.HasLength()
Base.length(acm::AbstractComMatrix) = size(acm, 1)

Base.propertynames(acm::AbstractComMatrix) = sitenames(acm)
Base.getproperty(acm::AbstractComMatrix, nm::Symbol) = acm[:, findfirst(==(nm), sitenames(acm))]

Tables.rows(acm::AbstractComMatrix) = acm
Tables.columns(acm::AbstractComMatrix) = eachcol(occurrences(acm))

Tables.schema(acm::AbstractComMatrix) = Tables.Schema(sitenames(acm), eltype.(eachcol(occurrences(acm))))


function ComMatrix(x)
    Tables.istable(x) || throw(ArgumentError("input is not a table"))
    rows = Tables.rows(x)
    sch = Tables.schema(rows)
    names = sch.names
    types = sch.types
    # custom constructor that creates an "empty" ComMatrix according to given column names & types
    # note that the "unknown" schema case should be considered, i.e. when `Tables.schema(x) === nothing`
    mytbl = ComMatrix(names, types)
    for row in rows
        # a convenience function provided in Tables.jl for "unrolling" access to each column/property of a `Row`
        # it works by applying a provided function to each value; see `?Tables.eachcolumn` for more details
        Tables.eachcolumn(sch, row) do val, columnindex::Int, columnname::Symbol
            push!(mytbl[columnindex], val)
        end
    end
    return mytbl
end

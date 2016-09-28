using SpatialEcology
using Base.Test

# write your own tests here
@test 1 == 1

# testing the constructor
function tes{T <: Union{Int, Bool}}(x::AbstractMatrix{T})
    println("success")
end

using SpatialEcology
using DataFrames
using CSV
using SparseArrays
using Test
import Random

Random.seed!(1337)

include("ComMatrix_tests.jl")
include("Assemblage_tests.jl")

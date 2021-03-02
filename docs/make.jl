using Documenter, SpatialEcology

makedocs(
    sitename="Spatial Ecology"
)

deploydocs(
    repo = "github.com/EcoJulia/SpatialEcology.jl.git",
    push_preview = true
)
using Documenter, SpatialEcology

makedocs(
    sitename="Spatial Ecology in Julia"
)

deploydocs(
    repo = "github.com/EcoJulia/SpatialEcology.jl.git",
    push_preview = true
)
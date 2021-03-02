using Documenter, SpatialEcology

makedocs(
    sitename="Spatial Ecology in Julia",
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/EcoJulia/SpatialEcology.jl.git",
    push_preview = true
)
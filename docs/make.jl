using Documenter, SpatialEcology

makedocs(
    sitename = "Spatial Ecology in Julia",
    authors = "Michael Krabbe Borregaard",
    pages = [
        "Home" => "index.md",
        "Tutorial / Quick start" => "tutorial.md",
        "Manual" => Any[
            "Creating SpatialEcology objects" => "man/constructors.md",
            "Working with SpatialEcology objects" => "man/getters.md",
            "Getting and setting data" => "man/data.md",
            "Subsetting" => "man/subsetting.md",
            "Randomizations" => "man/randomization.md",
        ],
        "Example analyses" => Any[
            "Node-based analysis of species distributions" => "examples/nodebased.md"
        ],
        "Library" => "lib/public.md"
    ]
)
 
deploydocs(
    repo = "github.com/EcoJulia/SpatialEcology.jl.git",
    push_preview = true
)
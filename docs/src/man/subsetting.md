# Groups and subsets

```@contents
Pages = ["subsetting.md"]
```

One of the most powerful ideas in SpatialEcology is that it lets you create views into all objects
(most importantly `ComMatrix` and `Assemblage`) based on a subset of species or sites. 
The object will drop unused species or sites.

Let's say for instance we want to calculate the average range size for each latitudinal band 
for the dataset of European amphibians.

First, we load the data:
```@example subset
using SpatialEcology, Plots, CSV, DataFrames, Statistics
ENV["GKSwstype"]="nul"
amphdata = CSV.read(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv"), DataFrame)
amph = Assemblage(amphdata[!, 4:end],amphdata[!, 1:3], sitecolumns = false);
```

And let's add the rangesizes of each species to the dataset
```@example subset
addtraits!(amph, occupancy(amph), :rangesize)
```

Then let's get all unique latitudes 
```@example subset
latitudes = unique(coordinates(amph)[:, 2])
```

We can use a simple to loop over all the latitudes, generate a relevant subset and calculate the mean rangesize
```@example subset
latitude_range = zeros(size(latitudes))
for (i, lat) in enumerate(latitudes)
    sites = findall(==(lat), coordinates(amph)[:,2])
    subset = view(amph, sites = sites)
    latitude_range[i] = mean(subset[:rangesize])
end
scatter(latitudes, latitude_range, xlab = "Latitude", ylab = "Mean range size")
```

Subsetting and sampling over a factor is common enough that there is a specialized syntax 
for this, `groupspecies` and `groupsites`. 
All of the above can be expressed by grouping the assemblage over the second coordinate (latitude):
```@example subset
latitudinal_assemblages = groupsites(amph, coordinates(amph)[:,2], dropspecies = true)
latitude_range = [mean(lat[:rangesize]) for lat in latitudinal_assemblages]
```

You can also use subsetting to plot a single species:
```@example subset
spec = view(amph, species = ["_Bufo_bufo"])
plot(spec, title = "Common Toad", showempty = true, c = cgrad([:grey, :red], categorical = true))
```
#Todo make this work without wrapping `sp`





Note that `view` does not copy the data. To create a new, independent `Assemblage` from a view, use `copy`.

## Indexing a ComMatrix

A `ComMatrix` is a [DimensionalData](https://rafaqz.github.io/DimensionalData.jl/) array with
the dimensions `Species` and `Sites`, so it indexes like any matrix, and species and site names are
selected with `At`:
```julia
com = commatrix(amph)
com[1, 3]                           # the value for the first species in the third site
com[At("_Bufo_bufo"), :]            # that species in all sites, as a vector named by site
com[1:10, [1, 3]]                   # a ComMatrix with the first ten species in two sites
com[Species = At(["_Bufo_bufo", "_Rana_temporaria"])]
```
As for an `Array`, a single species or site drops that dimension (the dropped name is kept in
the result's `refdims`); use a range, such as `com[1:1, :]`, to keep a one-row `ComMatrix`.
`@view com[1:10, :]` and `view(com, species = names)` give a `SubComMatrix` that shares the data.

You can also write into a `ComMatrix`, e.g. `com[1, 3] = true`, directly or through a view.
Writing a zero removes the record.


## Index

```@index
Pages = ["subsetting.md"]
```

## API
```@docs
groupspecies
groupsites
coarsen
```

### Utilities
```@docs
dispersionfield
pairwise
asquantiles
asquantiles!
```
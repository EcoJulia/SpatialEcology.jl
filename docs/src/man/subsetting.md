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
```@example
using SpatialEcology, Plots, CSV, DataFrames, Statistics
amphdata = CSV.read(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv"), DataFrame)
amph = Assemblage(amphdata[!, 4:end],amphdata[!, 1:3], sitecolumns = false);
```

And let's add the rangesizes of each species to the dataset
```@example
addtraits!(amph, occupancy(amph), :rangesize)
```

Then let's get all unique latitudes 
```@example
latitudes = unique(coordinates[:, 2])
```

We can use a simple to loop over all the latitudes, generate a relevant subset and calculate the mean rangesize
```@example
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
```@example
latitudinal_assemblages = groupsites(amph, coordinates(amph)[:,2], dropspecies = true)
latitude_range = [mean(lat[:rangesize]) for lat in latitudinal_assemblages]
```

You can also use subsetting to plot a single species:
```@example
spec = view(amph, species = ["_Bufo_bufo"])
plot(spec, title = "Common Toad", showempty = true, c = cgrad([:grey, :red], categorical = true))
```
#Todo make this work without wrapping `sp`





Note that `getindex` (`[]`) will create a view by default - to create a new `Assemblage` object you can use `copy`.


## Index

```@index
Pages = ["subsetting.md"]
```

## API
```@docs
groupspecies
groupsites
aggregate
```

### Utilities
```@docs
dispersionfield
pairwise
asquantiles
asquantiles!
```
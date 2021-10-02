# Getting started:

The aim of SpatialEcology is to make it easy to input species occurrence data, along with species traits and site attributes, 
and do most analyses intuitively and simply. For this example we'll use a dataset of the distributions of all amphibians
in Europe, based on a (previous) rasterization of IUCN's amphibian shapefile data onto a 0.5x0.5 lat-long grid (provided 
with the package).

First, let's load the relevant libraries
```@example tutorial
using SpatialEcology, Plots, CSV, DataFrames, Statistics
ENV["GKSwstype"]="nul"
```

We read in the species occurrence data from a DataFrame.
The object constructors take a wide range of input data types, a typical being a presence-absence matrix 
as a DataFrame along with the spatial coordinates of sites as a 3-column DataFrame.
In the example data, the site coordinates are simply the first three columns
```@example tutorial
amphdata = CSV.read(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv"), DataFrame);
amphdata[1:3,1:6]
```

We split the first three columns off the DataFrame and create an `Assemblage` object
The `sitecolumns` keyword tells SpatialEcology that the input DataFrame has sites as rows (and species as columns)
```@example tutorial
amph = Assemblage(amphdata[!, 4:end],amphdata[!, 1:3], sitecolumns = false)
```

Plotting recipes are defined for `Assemblage` objects - here we can simply plot the entire object to get a richness map
```@example tutorial
plot(amph)
```

## Accessing and adding data:
Access functions summarize the data, such as `occupancy`, `richness`, `nsites`, `nspecies`
```@example tutorial
mean(occupancy(amph))
```
Add DataFrames or Vectors of data to the assemblage, DataFrames are automatically aligned keeping everything together. Access data by column name.
```@example tutorial
addtraits!(amph, occupancy(amph), :rangesize)
histogram(amph[:rangesize], grid = false, legend = false)
```

## Easy subsetting and quick views:
It is easy to do analyses by taking subsets of the objects and analyzing on them - i.e.
the split-apply-combine strategy. Let's say we want to get the mean range of all species 
occurring in each grid cell. We can `map` over all sites in the assemblage, and for each
site use `occurring` to extract the indexes of the species that occur in that site. Finally
we can use this to index into the `rangesize` column we added to the dataframe above, and 
take the mean. A plot recipe allows us to create a heatmap for any numeric vector of values
that match the sites in the `Assemblage` object.
```@example tutorial
meanrange = map(site->mean(amph[:rangesize][occurring(amph,site)]), 1:nsites(amph))
plot(meanrange, amph, color = :fire)
```

The subsetting is fast as it uses `view`s by default and allocate very little. You
can also use `view` explicitly to get a handle on a subset of an `Assemblage`, subsetting
by either a vector of species (to get a smaller taxonomic group) or a vector of sites (to
get a smaller geographic area). Here we get the names of all species in the genus `Triturus`
and use this to create a view that acts as a smaller `Assemblage` for this species alone.
```@example tutorial
triturus = view(amph, species = contains.(speciesnames(amph), "Triturus"))
```
We can then use this dataset for further analyses - here getting the latitudinal range for
the genus in Europe:
```@example tutorial
extrema(coordinates(triturus)[:,1])
```

## Aggregation and other operations
For gridded Assemblages, you can also do some simple geographic operations on the object,
such as aggregrating the grid to a coarser grain size. Here, we lump them to a 2 times
coarser grain
```@example tutorial
amp2 = aggregate(amph, 2)
default(color = cgrad(:Spectral, rev = true))
plot(plot(amph), plot(amp2))
```

There are currently a few extra less generic operations defined on Assemblages, that
may at some point be moved to a different module. For example, you can generate the 
"dispersion field" of a given site using the `dispersionfield` function. A dispersion
field of a site is a map of the total richness of all species occurring in a given site 
[Graves & Rahbek 2005 PNAS]("http://macroecointern.dk/pdf-reprints/Graves_and_Rahbek_PNAS_2005.pdf").
This illustrates the compositional similarity of the site to it's surrounding sites. Here
we plot the dispersion field for site number 50. 
```@example tutorial
plot(dispersionfield(amph, 50), amph, c = :rainbow)
```
This functionality of the SpatialEcology package was used for the analyses in [Borregaard, Graves & Rahbek, 2020]("http://macroecointern.dk/pdf-reprints/Borregaard_2020_Nature.pdf").


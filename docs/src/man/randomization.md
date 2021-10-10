# Randomization

```@contents
Pages = ["randomization.md"]
```

The methods for matrix randomizations defined in 
[RandomBooleanMatrices](https://github.com/ecojulia/RandomBooleanMatrices.jl)
are defined for the SpatialEcology type `ComMatrix` and `Assemblage`. These 
methods randomize the cells of a presence-absence matrix while keeping row and 
column sums constant - a very widely used method for randomizing ecological
communities for null model analysis (see e.g. 
[Miklós, I. and Podani, J. (2004), RANDOMIZATION OF PRESENCE–ABSENCE MATRICES: COMMENTS AND NEW ALGORITHMS. Ecology, 85: 86-92](https://doi.org/10.1890/03-0101).

To simulate random communities, you define a `matrixrandomizer` object on your 
type - you can then simulate from it using `rand` and `rand!`.
For example, using the `amph` object from the tutorial
```@example tutorial
using Random
rmg = matrixrandomizer(amph)
newcomm = rand!(rmg)
````

`newcomm` is then a new `Assemblage` object with the same properties as `amph`, but with all co-occurrences randomized. 


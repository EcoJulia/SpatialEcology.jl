# SpatialEcology

[![d_stable](https://img.shields.io/badge/Doc-stable-green?style=flat-square)](https://ecojulia.github.io/SpatialEcology.jl/stable/)
[![d_latest](https://img.shields.io/badge/Doc-latest-blue?style=flat-square)](https://ecojulia.github.io/SpatialEcology.jl/dev/)

![version](https://img.shields.io/github/v/tag/EcoJulia/SpatialEcology.jl?sort=semver&style=flat-square)
[![CI](https://github.com/EcoJulia/SpatialEcology.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/EcoJulia/SpatialEcology.jl/actions/workflows/CI.yml)
![Doc](https://img.shields.io/github/workflow/status/EcoJulia/SpatialEcology.jl/Documentation?label=Doc&style=flat-square)
[![codecov](https://codecov.io/gh/EcoJulia/SpatialEcology.jl/graph/badge.svg?token=DeSFZuHa99)](https://codecov.io/gh/EcoJulia/SpatialEcology.jl)

## Primary author: Michael Krabbe Borregaard (@mkborregaard)

A package for community- and macro-ecological analysis in julia.
This package offers a set of base types for interoperability in spatial ecology. The idea is to provide a powerful framework for expressing a great variety of analyses in a flexible manner. It presently holds types for presence-absence matrices, site data and species traits, and will be included with phylogenies and ecological interaction networks. SpatialEcology takes care of aligning all data for analysis.

The emphasis is on fast, flexible code operating mainly with views into the larger dataset. It currently holds fast, specialized code for operations on views into sparse matrices (such as presence-absence matrices). This allows analyses to be done in a split-apply-combine framework.

The package originated as a port of the R package `nodiv`, available from CRAN.

- Types:
  - Assemblage (holds presence-absence information along with information on traits and sites)
  - ComMatrix (presence-absence matrix)
  - SpatialData (Grid or Point data with site information)

## Relevant other packages

This package is part of the [EcoJulia](https://ecojulia.org) organisation, which aims to bring together a coherent set of packages for ecological data analysis.For other relevant packages check the [BioJulia](https://biojulia.net) organisation focusing on molecular biology, and the [JuliaGeo](https://juliageo.org/) organisation focusing on geographical data analysis. A long-term goal of the EcoJulia organisation is to interface as seamlessly as possible with these projects to create an integrated data analysis framework for julia.

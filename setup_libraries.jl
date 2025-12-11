using Pkg

# standard packages (strings)
Pkg.add([
    "DelimitedFiles",
    "Plots",
    "Measures",
    "StatsBase",
    "DataFrames",
    "StatsPlots",
    "OrdinaryDiffEq",
    "LsqFit",
    "NonlinearSolve",
    "Statistics",
])

# add the GitHub repo as a package
Pkg.add(PackageSpec(url="https://github.com/ComputationalMechanobiology/CellAdhesion.jl.git"))

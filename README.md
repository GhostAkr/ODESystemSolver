# ODE system solver

This repository contains a set of useful algorithms which can be used to solve a system of 
ODEs.

## Repository structure

- `RungeKutta.jl` contains Runge - Kutta algorithm. Order and number of total stages can be set by user. Method's coefficients are not calculated and should be passed to the function as arguments.

- `DormandPrince.jl` contains Dormand - Prince algorithm.

- `Task.jl` contains problem examples (system of ODEs).

- `Main.jl` can be considered as an entry-point from which you can run tests from `Task.jl`.

- `Charts.nb` contains charts in Wolfram Mathematica analysing obtained results.

## Example

As a test example we consider following problem (defined in `Task.jl`): [Task Problem](docs/example_problem.pdf).


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

As a test example we consider following system of ODE (defined in `Task.jl`):

$$
y'_1 = 2t \cdot y_1 \cdot \text{ln}(\text{max}(y_2, 10^{-3})),
$$

$$
y'_2 = -2t \cdot y_2 \cdot \text{ln}(\text{max}(y_1, 10^{-3})).
$$

Known analytical solution for this problem is:

$$
y_1(t) = \exp(\sin(t^2)),
$$

$$
y_2(t) = \exp(\cos(t^2)).
$$

Integration interval is $[0.1, 4.1]$. As initial conditions we consider analytical solutions
at $t_1$: $y_1(t_1)$ and $y_2(t_1)$.

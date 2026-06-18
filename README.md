# Backward Difference Method for Parabolic Partial Differential Equations

A Python implementation of an implicit finite-difference method for the one-dimensional heat equation. The project discretizes the spatial domain with centered differences, advances the solution with a backward difference in time, and solves the resulting linear system through a generalized Crout factorization routine.

## Mathematical problem

The code considers the parabolic partial differential equation

$$\frac{\partial u}{\partial t}(x,t)=\alpha^2\frac{\partial^2u}{\partial x^2}(x,t),\qquad 0<x<l,\quad t>0,$$

subject to homogeneous Dirichlet boundary conditions and an initial condition:

$$u(0,t)=u(l,t)=0,\qquad u(x,0)=f(x).$$

The spatial and temporal step sizes are

$$h=\frac{l}{m},\qquad k=\frac{T}{N},$$

with grid points

$$x_i=ih,\qquad t_j=jk.$$

## Backward-difference discretization

The time derivative is approximated implicitly at time level $t_j$:

$$\frac{\partial u}{\partial t}(x_i,t_j)\approx\frac{u(x_i,t_j)-u(x_i,t_{j-1})}{k}.$$

The second spatial derivative is approximated with a centered difference:

$$\frac{\partial^2u}{\partial x^2}(x_i,t_j)\approx\frac{u(x_{i+1},t_j)-2u(x_i,t_j)+u(x_{i-1},t_j)}{h^2}.$$

After introducing

$$\lambda=\frac{\alpha^2k}{h^2},$$

the numerical scheme becomes

$$(1+2\lambda)w_{i,j}-\lambda w_{i+1,j}-\lambda w_{i-1,j}=w_{i,j-1}.$$

At every time step, the interior solution therefore satisfies

$$A\mathbf{w}^{(j)}=\mathbf{w}^{(j-1)},$$

where the coefficient matrix can be written compactly as

$$A=\mathrm{tridiag}\!\left(-\lambda,\,1+2\lambda,\,-\lambda\right).$$

For $\lambda>0$, $A$ is symmetric, positive definite, and strictly diagonally dominant. The backward-difference method is unconditionally stable, and its truncation error is

$$O(k+h^2).$$

## Repository files

| File | Description |
|---|---|
| `backward_difference_method.py` | Builds the finite-difference system, advances the numerical solution, and exports an error table. |
| `crout_factorization_generalization.py` | Implements a generalized Crout solver for block-tridiagonal systems. |
| `Backward_Difference_Parabolic.tex` | Original Spanish mathematical report. |
| `Backward_Difference_Parabolic_EN.tex` | English translation of the report. |

## Numerical example

The included example solves

$$\frac{\partial u}{\partial t}-\frac{\partial^2u}{\partial x^2}=0,\qquad 0<x<1,$$

with

$$u(0,t)=u(1,t)=0,\qquad u(x,0)=\sin(\pi x).$$

The analytical solution is

$$u(x,t)=e^{-\pi^2t}\sin(\pi x).$$

The default parameters are:

| Parameter | Value |
|---|---:|
| Spatial endpoint `l` | 1 |
| Final time `T` | 0.5 |
| Diffusion coefficient `alpha` | 1 |
| Time subdivisions `N` | 50 |
| Spatial subdivisions `m` | 10 |
| Spatial step `h` | 0.1 |
| Time step `k` | 0.01 |

At the midpoint $x=0.5$ and final time $t=0.5$, the implementation produces the approximation

$$w_{5,50}\approx 0.00937818,$$

while the analytical value is

$$u(0.5,0.5)\approx 0.00719188.$$

The corresponding absolute error is approximately

$$2.18630\times10^{-3}.$$

## Requirements

- Python 3
- NumPy
- pandas

Install the dependencies with:

```bash
python -m pip install numpy pandas
```

## Usage

Place both Python modules in the same directory and run:

```bash
python backward_difference_method.py
```

The script prints the numerical state at each time level and creates:

```text
error_table.csv
```

The CSV file contains the spatial coordinate, numerical approximation, analytical solution, and absolute error at $t=0.5$.

## Implementation outline

1. Construct the uniform spatial and temporal grids.
2. Evaluate the initial condition at the interior spatial points.
3. Assemble the tridiagonal coefficient matrix.
4. Solve the implicit system at each time step.
5. Compare the numerical result with the analytical solution.
6. Export the pointwise errors to a CSV file.

## Current implementation detail

The time array contains $N+1$ values, and the current loop performs a linear solve for every value in that array. The function stores the state immediately before the last solve in `aux`; consequently, the example uses the returned `w_aux` value when constructing the error table at $t=T$.

## Author

**BSc. Julio Medina**  
Universidad de San Carlos de Guatemala  
School of Physical and Mathematical Sciences  
Master's Program in Physics

## References

1. Richard L. Burden and J. Douglas Faires, *Numerical Analysis*, 9th edition, Brooks/Cole, Cengage Learning.
2. Richard S. Varga, *Matrix Iterative Analysis*, 2nd edition, Springer.
3. Julio Medina, *Finite Difference Method for Elliptic Equations*.

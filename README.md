# VintageOptimalControl

[![Build Status](https://github.com/michaelfreiberger/VintageOptimalControl.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/michaelfreiberger/VintageOptimalControl.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a Julia-package which allows for the solution of age-structured/vintage-structured optimal control problems.

# Installation

The package is available in the Julia repository. The package can be installed by simply using the following two commands:

```julia-repl
julia> using Pkg
julia> Pkg.add("VintageOptimalControl")
```

# Documentation

For information on how to use the toolbox and some examples, please take a look into the [documentation](https://michaelfreiberger.github.io/VintageOptimalControl.jl/dev/). 

# Problem class

This package is designed to solve age-structured problems of the form

$$ 
\max_{u: [0,T]\to\mathbb{R}^{n_u},\; v:[0,T]\times[0,\omega]\to\mathbb{R}^{n_v}}\int_{0}^{T}\int_{0}^{\omega} \exp(-\rho t)\cdot J\Big(X(t),Y(t,a),Z(t),u(t),v(t,a),t,a\Big) da dt + \int_0^{\omega} S(X(T),Y(T,a),Z(T),T,a)da
$$

$$
s.t. \hspace{2cm}\dot{X}(t) = f(X(t),Z(t),u(t),t) ,\hspace{3cm} X(0) = X_0
$$

$$
\hspace{3.8cm} \Big(\frac{\partial}{\partial t}+\frac{\partial}{\partial a}\Big) Y(t,a) = g(X(t),Y(t,a),Z(t),u(t),v(t,a),t,a) ,\qquad Y(t,0) = y_0(X(t),Z(t),u(t),t)
$$

$$
\hspace{0.5cm} Z(t) = \int_0^{\omega} h(X(t),Y(t,a),u(t),v(t,a),t,a) da
$$

$$
\underline{u_i} \leq u_i(t) \leq \overline{u_i} \qquad\forall i = 1,\ldots,n_u \quad,\quad\forall t\in [0,T]
$$

$$
\hspace{2.6cm}\underline{v_i} \leq v_i(t,a) \leq \overline{v_i} \qquad\forall i = 1,\ldots,n_v \quad,\quad\forall t\in [0,T]\quad,\quad \forall s\in [0,\omega]
$$

A vintage-structured optimal control problem can be easily transformed into a vintage-structured optimal control problem considering the relationship

$$
    s = t-a
$$
between vintage $s$, time $t$, and age $a$.

## **Variable and function explanations**

The model consist of three different types of variables:

* concentrated state variables $X$ and control variables $u$ (depending on time $t$)
* distributed state variables $Y$ and control variables $v$ (depending on time $t$ and age $a$)
* aggregated state variables $Z$ (depending on time $t$)

### Concentrated variables

The concentrated variables consist of the $n_u$-dimensional control variable $u$ and the $n_X$-dimensional state variable $X$. The initial value of $X$ is given by $X_0$ and its dynamics are affected by $X$, $Z$, and $u$. Furthermore $X$ and $u$ also affect the initial values of the distributed variables $Y$ (see below), the dynamics of $Y$ and the aggregated variables $Z$.

### Distributed variables

Analogous to the concentrated variables, the distributed variables consist of the $n_v$-dimensional control variable $v$ and the $n_Y$-dimensional state variable $Y$. The (strictly speaking partial) differential equation can be solved along the characteristic lines with fixed vintage $t-a(=s)=const$ . While the dynamics of $Y$ are affected by the distributed controls,states as well as the concentrated states and controls and aggregated states. The initial values of $Y$ for age $0$ at each time $t$ are determined by the value of the concentrated state $Z$ at time $t$, the concentrated control $u$, and the aggregated variable $Z$.

### Aggregated variables

The aggregated variable $Z$ at each point in time is defined through the aggregation over all ages. The function $h()$ captures how the concentrated and distributed variables are aggregated.


## **Algorithmic principles**

The algorithm uses a gradient based approach following a sequence of steps:

1. Start with a given guess for the optimal solution of the control variables.
2. Calculate the corresponding profiles of the state variables and co-state variables (based on the maximum principle)
3. Calculate the gradient based on the Hamiltonian.
4. Adjust the guess for the control in the direction of the gradient.
5. Find the optimal adjustment step in the direction of the gradient.
6. Define new currently best solution.
7. Iterate from step 2 until no improvement can be found anymore.

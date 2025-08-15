# Mandelbrot Geometry and Stability in Population Growth Models

This repository contains the code and research materials for my IB Mathematics HL Extended Essay.  
The project investigates the relationship between **period-n bulbs of the Mandelbrot set** and the **stability regions in the bifurcation diagram of the logistic map**.

---

## Research Question

**How do the period-n bulbs of the Mandelbrot set correspond to stability regions in the bifurcation diagram of the logistic map?**

---

## Background

### Mandelbrot Set

The Mandelbrot set is defined by the recursive relation:

$$
z_{n+1} = z_n^2 + c, \quad z_0 = 0, \quad c \in \mathbb{C}
$$

A complex number $c$ belongs to the Mandelbrot set if the sequence $\{z_n\}$ remains bounded as $n \to \infty$.  
Each **bulb** of the Mandelbrot set corresponds to parameter regions where the iteration settles into a cycle of given period.

### Logistic Map

The logistic map models population growth and is given by:

$$
x_{n+1} = r x_n (1 - x_n), \quad x_0 \in (0,1), \quad r \in (0,4]
$$

where:
- $r$ = growth rate parameter  
- $x_n$ = normalized population at time step $n$

The bifurcation diagram of the logistic map shows **stability windows** — intervals of $r$ where the system converges to a cycle of fixed period.

---

## Methodology

### Analytical Approach

For the logistic map:

Candidate cycle points of period $k$ are solutions of

$$
f_r ^{(k)}(x) = x
$$

where $f_r(x) = r x (1-x)$ and $f_r^{(k)}$ denotes the $k$-th iterate.

Stability is determined by:

$$
\left| (f_r^{(k)})'(x) \right| = \prod_{j=0}^{k-1} | f_r'(x_j) | < 1
$$

For the Mandelbrot set:

Period-$n$ bulbs correspond to parameters $c$ such that:

$$
f_c^{(n)}(0) = 0, \quad f_c^{(m)}(0) \neq 0 \quad \forall m < n
$$

with $f_c(z) = z^2 + c$.  
Stability of these cycles is confirmed when:

$$
\left| (f_c^{(n)})'(z^*) \right| < 1
$$

### Computational Approach

- **Bifurcation diagram**:
- Sweep $r \in [2.8, 4]$ on a fine grid.
- Iterate the logistic map, discard burn-in, and plot subsequent values.
- Detect cycle periods numerically.

- **Mandelbrot set**:
- Generate via escape-time algorithm with $R = 2$.
- Identify bulbs and verify period-$n$ dynamics.

- **Comparison**:
- Match stability windows (logistic map) with bulbs (Mandelbrot set).
- Pair diagrams side-by-side to illustrate correspondence.

---

## Results (Preview)

- Period-1 region ($1 < r < 3$) in logistic map ↔ Main cardioid of Mandelbrot set.  
- Period-2 stability window in logistic map ↔ Period-2 bulb in Mandelbrot set.  
- Higher-order bulbs (period-3, period-4, etc.) ↔ Corresponding stability regions in bifurcation diagram.

---

## Repository Contents

- `src/` — Python scripts for:
  - Generating Mandelbrot set visualisations
  - Computing and plotting logistic map bifurcation diagrams
- `notebooks/` — Jupyter notebooks with computational experiments
- `EE.pdf` — Full Extended Essay writeup


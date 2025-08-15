# Mandelbrot-LogisticMap-EE

Code and visualisations for my IB Mathematics AA HL Extended Essay.

**Research Question:** *How do the period-n bulbs of the Mandelbrot set correspond to stability regions in the bifurcation diagram of the logistic map?*

---

## Overview

This repo contains the analysis and code used to study a period-preserving link between:

- **The Mandelbrot set** (quadratic family)

  $$
  z_{n+1}=z_n^2+c
  $$

  where the interior “bulbs” correspond to attracting cycles of a given period.

- **The Logistic map** (population model)

  $$
  x_{n+1}=rx_n(1-x_n)
  $$

  whose bifurcation diagram shows windows of stable period-$k$ behaviour interleaved with chaos.

The goal is to **match period-$n$ Mandelbrot bulbs with period-$k$ stability windows of the logistic map** and show they align in period (i.e., $n=k$) via analytic conditions and computational experiments.

---

## Methodology

### Analytical

- **Logistic map cycles and stability**

  $$
  f_r(x)=r\,x(1-x),\qquad f_r^{\,k}(x)=x,\qquad
  \left|\left(f_r^{\,k}\right)'(x)\right|=\prod_{j=0}^{k-1}\left|f_r'(x_j)\right|
  =\prod_{j=0}^{k-1}\left|r(1-2x_j)\right|<1
  $$

  Bifurcation boundaries occur where $\left|\left(f_r^{\,k}\right)'(x)\right|=1$.

- **Mandelbrot bulbs (quadratic family)**

  $$
  f_c(z)=z^2+c,\qquad
  f_c^{\,n}(0)=0 \ \text{ and }\ f_c^{\,m}(0)\neq 0 \ \ (\forall\, m<n)
  $$

  Centers of period-$n$ bulbs are superattracting; in the interior the multiplier of the $n$-cycle satisfies $| (f_c^{\,n})'(z^\*) |<1$.

- **Link:** match the integer period $n$ (bulb) with $k$ (logistic window) and verify via the conditions above.

### Computational

- **Bifurcation diagram**
  1. Sweep $r$ (e.g., $[2.8,4]$), iterate from $x_0\in(0,1)$.
  2. Discard burn-in, plot subsequent iterates.
  3. Detect periods via clustering/repeats to label stability windows.

- **Mandelbrot set**
  1. Escape-time algorithm on a grid of $c$ with $z_0=0$, radius $R=2$, max iters $N$.
  2. Locate/label period-$n$ bulbs; verify multipliers $<1$.

- **Comparison**
  - Place labelled bifurcation windows next to labelled Mandelbrot bulbs; tabulate $(k,\ r\text{-interval}) \leftrightarrow (n,\ c\text{-region})$.

---

## Features

- High-res **bifurcation diagrams** with optional period labelling.  
- **Mandelbrot renderer** with zoom and optional bulb highlighting.  
- Side-by-side comparison utilities and figure exporters for the EE.

---

## Quick Start

```bash
git clone https://github.com/matanou/Mandelbrot-LogisticMap-EE.git
cd Mandelbrot-LogisticMap-EE

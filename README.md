# Mandelbrot-LogisticMap-EE

Code and visualisations for my IB Mathematics AA HL Extended Essay.

**Research Question:** *How do the period-n bulbs of the Mandelbrot set correspond to stability regions in the bifurcation diagram of the logistic map?*

> **Note on equations:** GitHub READMEs don’t render LaTeX natively.  
> All formulas below are embedded as SVG images with a white background so they display clearly in both dark and light themes.

---

## Overview

This repo contains the analysis and code used to study a period-preserving link between:

- **The Mandelbrot set** (quadratic family)

  ![z recurrence](https://latex.codecogs.com/svg.latex?\bg_white%20z_%7Bn%2B1%7D%20%3D%20z_n%5E2%20%2B%20c)

  where the interior “bulbs” correspond to attracting cycles of a given period.

- **The Logistic map** (population model)

  ![logistic map](https://latex.codecogs.com/svg.latex?\bg_white%20x_%7Bn%2B1%7D%20%3D%20r\,x_n(1-x_n))

  whose bifurcation diagram shows windows of stable period-\(k\) behaviour interleaved with chaos.

The goal is to **match period-\(n\) Mandelbrot bulbs with period-\(k\) stability windows of the logistic map** and show they align in period (i.e., \(n=k\)) via analytic conditions and computational experiments.

---

## Methodology

### Analytical

- **Logistic map cycles and stability**

  ![cycle and multiplier](https://latex.codecogs.com/svg.latex?\bg_white%20f_r(x)%3D%20r\,x(1-x),\quad%20f_r%5E{k}(x)%3D%20x,\quad%20\left%7C\left(f_r%5E{k}\right)'(x)\right%7C%20%3D%20\prod_%7Bj%3D0%7D%5E%7Bk-1%7D\left%7C%20r(1-2x_j)\right%7C%20%3C%201)

  Bifurcation boundaries occur where

  ![boundary](https://latex.codecogs.com/svg.latex?\bg_white%20\left%7C%20\left(f_r%5E{k}\right)'(x)%20\right%7C%20%3D%201.)

- **Mandelbrot bulbs (quadratic family)**

  ![bulb condition](https://latex.codecogs.com/svg.latex?\bg_white%20f_c(z)%3D%20z%5E2%2Bc,\quad%20f_c%5E{n}(0)%3D0%20\text{%20and%20}%20f_c%5E{m}(0)\neq0\%20(\forall\,m%3Cn))

  Centers are superattracting; interior points have multiplier \(|(f_c^{\,n})'(z^\*)|<1\).

- **Link:** match period \(n\) (bulb) with \(k\) (logistic window) and verify via the conditions above.

### Computational

- **Bifurcation diagram**
  1. Sweep \(r\) (e.g., \([2.8,4]\)), iterate from \(x_0\in(0,1)\).
  2. Discard burn-in, plot subsequent iterates.
  3. Detect periods via clustering/repeats to label stability windows.

- **Mandelbrot set**
  1. Escape-time algorithm on a grid of \(c\) with \(z_0=0\), radius \(R=2\), max iters \(N\).
  2. Locate/label period-\(n\) bulbs; verify multipliers \(<1\).

- **Comparison**
  - Place labelled bifurcation windows next to labelled Mandelbrot bulbs; tabulate \((k,\ r\text{-interval}) \leftrightarrow (n,\ c\text{-region})\).

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

# Figures

All figures are self-generated from code in `src/`. Each figure is reproducible via the exact commands below.
Image files are saved at high DPI and without axes so they can be dropped into LaTeX with a white background.

> **Provenance note**: These images are "Author’s own visualisations" generated for the Extended Essay.
Repro steps and parameters are recorded here for transparency and IB academic integrity.

---

## Figure Index

| Label in EE | Filename                      | Source script                | How to reproduce (example command)                                             |
|-------------|-------------------------------|------------------------------|--------------------------------------------------------------------------------|
| Figure 1    | `mandelbrot.png`              | `src/mandelbrot_vis.py`      | see command below                                                              |
| Figure 2    | `bifurcation.png`             | `src/bifurcation_of_logistic_map.py`| see command below                                                              |

---

## 1) Mandelbrot set visualisation (Figure 1)

**Purpose**: Visualise membership in the Mandelbrot set by escape-time colouring of points within radius 2 of the origin, consistent with the EE text.  
**Default region**: `xmin=-2.5, xmax=1.0, ymin=-1.5, ymax=1.5` (full set framing)  
**Output**: `figures/mandelbrot.png`

```bash
python src/mandelbrot_vis.py \
  --xmin -2.5 --xmax 1.0 --ymin -1.5 --ymax 1.5 \
  --width 2400 --height 1800 --max-iter 500 \
  --cmap magma \
  --out figures/mandelbrot.png

python src/logistic_bifurcation.py \
  --r-min 2.5 --r-max 4.0 --r-steps 5000 \
  --x0 0.5 \
  --n-transient 1000 --n-keep 200 \
  --out figures/bifurcation.png

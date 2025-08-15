
---

## 4) `src/mandelbrot_vis.py`

```python
#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Optional speed-up via numba if present
try:
    from numba import njit, prange
    NUMBA = True
except Exception:
    NUMBA = False

def _mandelbrot_numpy(xmin, xmax, ymin, ymax, width, height, max_iter):
    xs = np.linspace(xmin, xmax, width, dtype=np.float64)
    ys = np.linspace(ymin, ymax, height, dtype=np.float64)
    image = np.zeros((height, width), dtype=np.uint16)
    for iy, y0 in enumerate(ys):
        for ix, x0 in enumerate(xs):
            x, y, it = 0.0, 0.0, 0
            while x*x + y*y <= 4.0 and it < max_iter:
                x, y = x*x - y*y + x0, 2.0*x*y + y0
                it += 1
            image[iy, ix] = it
    return image

if NUMBA:
    @njit(parallel=True, fastmath=True)
    def _mandelbrot_numba(xmin, xmax, ymin, ymax, width, height, max_iter):
        image = np.zeros((height, width), dtype=np.uint16)
        for iy in prange(height):
            y0 = ymin + (ymax - ymin) * iy / (height - 1)
            for ix in range(width):
                x0 = xmin + (xmax - xmin) * ix / (width - 1)
                x = 0.0
                y = 0.0
                it = 0
                while x*x + y*y <= 4.0 and it < max_iter:
                    x, y = x*x - y*y + x0, 2.0*x*y + y0
                    it += 1
                image[iy, ix] = it
        return image

def mandelbrot(xmin, xmax, ymin, ymax, width, height, max_iter, use_numba=True):
    if use_numba and NUMBA:
        return _mandelbrot_numba(xmin, xmax, ymin, ymax, width, height, max_iter)
    return _mandelbrot_numpy(xmin, xmax, ymin, ymax, width, height, max_iter)

def render(image, out, cmap="magma", dpi=200, transparent=False):
    plt.figure(figsize=(image.shape[1]/dpi, image.shape[0]/dpi), dpi=dpi)
    plt.imshow(image, extent=[0,1,0,1], cmap=cmap, origin="lower")
    plt.axis("off")
    plt.tight_layout(pad=0)
    plt.savefig(out, bbox_inches="tight", pad_inches=0, transparent=transparent)
    plt.close()

def parse_args():
    p = argparse.ArgumentParser(description="Generate a Mandelbrot render (PNG).")
    p.add_argument("--xmin", type=float, default=-2.5)
    p.add_argument("--xmax", type=float, default=1.0)
    p.add_argument("--ymin", type=float, default=-1.5)
    p.add_argument("--ymax", type=float, default=1.5)
    p.add_argument("--width", type=int, default=1200)
    p.add_argument("--height", type=int, default=900)
    p.add_argument("--max-iter", type=int, default=300)
    p.add_argument("--cmap", type=str, default="magma")
    p.add_argument("--no-numba", action="store_true", help="Force pure NumPy (debug/compat).")
    p.add_argument("--out", type=str, default="figures/mandelbrot.png")
    return p.parse_args()

def main():
    args = parse_args()
    img = mandelbrot(
        args.xmin, args.xmax, args.ymin, args.ymax,
        args.width, args.height, args.max_iter,
        use_numba=not args.no_numba
    )
    render(img, args.out, cmap=args.cmap)
    print(f"Saved Mandelbrot image -> {args.out}")

if __name__ == "__main__":
    main()

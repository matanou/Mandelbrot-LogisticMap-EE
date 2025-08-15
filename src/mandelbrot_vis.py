import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import io, time, threading

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector

# Try to use numba for speed
try:
    from numba import njit, prange
    NUMBA = True
except ImportError:
    NUMBA = False

# High precision library
from mpmath import mp


# ---------------------- Mandelbrot kernels ----------------------
if NUMBA:
    @njit(parallel=True, fastmath=True)
    def mandelbrot_calc(xmin, xmax, ymin, ymax, width, height, max_iter):
        image = np.zeros((height, width), dtype=np.uint16)
        for i in prange(height):
            y0 = ymin + (ymax - ymin) * i / height
            for j in range(width):
                x0 = xmin + (xmax - xmin) * j / width
                x, y = 0.0, 0.0
                it = 0
                while x*x + y*y <= 4.0 and it < max_iter:
                    x, y = x*x - y*y + x0, 2.0*x*y + y0
                    it += 1
                image[i, j] = it
        return image
else:
    def mandelbrot_calc(xmin, xmax, ymin, ymax, width, height, max_iter):
        image = np.zeros((height, width), dtype=np.uint16)
        for i in range(height):
            y0 = ymin + (ymax - ymin) * i / height
            for j in range(width):
                x0 = xmin + (xmax - xmin) * j / width
                x, y = 0.0, 0.0
                it = 0
                while x*x + y*y <= 4.0 and it < max_iter:
                    x, y = x*x - y*y + x0, 2.0*x*y + y0
                    it += 1
                image[i, j] = it
        return image


def mandelbrot_highprec(xmin, xmax, ymin, ymax, width, height, max_iter, dps):
    """Arbitrary-precision fallback for ultra-deep zooms."""
    mp.dps = int(dps)
    image = np.zeros((height, width), dtype=np.uint16)
    xmin_m = mp.mpf(xmin); xmax_m = mp.mpf(xmax)
    ymin_m = mp.mpf(ymin); ymax_m = mp.mpf(ymax)
    dx = (xmax_m - xmin_m) / width
    dy = (ymax_m - ymin_m) / height
    four = mp.mpf(4)
    for i in range(height):
        y0 = ymin_m + dy * i
        for j in range(width):
            x0 = xmin_m + dx * j
            x = mp.mpf(0); y = mp.mpf(0)
            it = 0
            while x*x + y*y <= four and it < max_iter:
                x, y = x*x - y*y + x0, 2*x*y + y0
                it += 1
            image[i, j] = it
    return image


# ---------------------- App ----------------------
class MandelbrotApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Mandelbrot Explorer (Tkinter)")
        self.geometry("1200x800")
        self.minsize(1000, 700)

        # ----- Variables -----
        # Style
        self.dark_mode = tk.BooleanVar(value=True)
        self.enable_color = tk.BooleanVar(value=True)
        self.color_map = tk.StringVar(value="turbo")

        # Render params
        self.width = tk.IntVar(value=800)
        self.height = tk.IntVar(value=600)
        self.max_iter = tk.IntVar(value=300)

        # Bounds & baseline
        self.xmin = tk.DoubleVar(value=-2.0)
        self.xmax = tk.DoubleVar(value=1.0)
        self.ymin = tk.DoubleVar(value=-1.2)
        self.ymax = tk.DoubleVar(value=1.2)
        self._defaults = (-2.0, 1.0, -1.2, 1.2)
        self._start_span = (self._defaults[1] - self._defaults[0],
                            self._defaults[3] - self._defaults[2])

        # Zoom / selector
        self.drag_zoom = tk.BooleanVar(value=True)
        self.zoom_factor_ctrl = tk.DoubleVar(value=2.0)  # factor per drag
        self.zoom_level = tk.DoubleVar(value=1.0)        # readout (× from start)

        # Iteration scaling
        self.auto_iter = tk.BooleanVar(value=False)
        self.iter_power = tk.DoubleVar(value=1.0)    # exponential power
        self._base_iter_snapshot = self.max_iter.get()  # captured when Generate pressed

        # High precision
        self.high_prec = tk.BooleanVar(value=False)
        self.decimal_places = tk.IntVar(value=80)

        # ----- Layout -----
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        self.sidebar = ttk.Frame(self, padding=10)
        self.sidebar.grid(row=0, column=0, sticky="nsw")
        self.main = ttk.Frame(self, padding=(10, 10, 10, 10))
        self.main.grid(row=0, column=1, sticky="nsew")
        self.main.columnconfigure(0, weight=1)
        self.main.rowconfigure(0, weight=1)

        # Figure
        self.fig = Figure(figsize=(7, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.main)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        # Status
        self.status = ttk.Label(self.main, text="Ready.", anchor="w")
        self.status.grid(row=1, column=0, sticky="ew", pady=(8, 0))

        # Sidebar controls
        self._build_sidebar()

        # State
        self.image_data = None
        self.selector = None  # RectangleSelector instance

        # First render
        self._apply_dark_mode()
        self.after(50, self._generate_async)

    # -------- Sidebar --------
    def _build_sidebar(self):
        self.sidebar.columnconfigure(0, weight=0)
        self.sidebar.columnconfigure(1, weight=1)

        # Parameters
        lf = ttk.LabelFrame(self.sidebar, text="Mandelbrot Parameters", padding=10)
        lf.grid(row=0, column=0, columnspan=2, sticky="new", pady=(0, 10))

        def row(r, label, widget):
            ttk.Label(lf, text=label).grid(row=r, column=0, sticky="w", padx=(0, 6))
            widget.grid(row=r, column=1, sticky="ew")

        row(0, "Width", ttk.Spinbox(lf, from_=100, to=8000, textvariable=self.width, increment=100))
        row(1, "Height", ttk.Spinbox(lf, from_=100, to=8000, textvariable=self.height, increment=100))
        row(2, "Max Iter", ttk.Spinbox(lf, from_=10, to=20000, textvariable=self.max_iter, increment=10))
        row(3, "xmin", ttk.Entry(lf, textvariable=self.xmin))
        row(4, "xmax", ttk.Entry(lf, textvariable=self.xmax))
        row(5, "ymin", ttk.Entry(lf, textvariable=self.ymin))
        row(6, "ymax", ttk.Entry(lf, textvariable=self.ymax))

        # Style / Actions
        sf = ttk.LabelFrame(self.sidebar, text="Style / Actions", padding=10)
        sf.grid(row=1, column=0, columnspan=2, sticky="new", pady=(0, 10))
        for c in (0, 1): sf.columnconfigure(c, weight=1 if c == 1 else 0)

        ttk.Checkbutton(sf, text="Dark mode", variable=self.dark_mode,
                        command=self._apply_dark_mode).grid(row=0, column=0, columnspan=2, sticky="w")

        ttk.Checkbutton(sf, text="Enable color", variable=self.enable_color,
                        command=self._redraw_last).grid(row=1, column=0, columnspan=2, sticky="w")
        ttk.Label(sf, text="Color map").grid(row=2, column=0, sticky="w", padx=(0, 6))
        cmap_box = ttk.Combobox(sf, textvariable=self.color_map, state="readonly",
                                values=["turbo", "plasma", "inferno", "viridis", "twilight_shifted", "magma", "cividis"])
        cmap_box.grid(row=2, column=1, sticky="ew")
        cmap_box.bind("<<ComboboxSelected>>", lambda e: self._redraw_last())

        # Drag-zoom
        ttk.Checkbutton(sf, text="Drag-zoom (select → center zoom)",
                        variable=self.drag_zoom, command=self._toggle_selector)\
            .grid(row=3, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Label(sf, text="Zoom factor").grid(row=4, column=0, sticky="w", padx=(0, 6))
        ttk.Spinbox(sf, from_=1.2, to=20.0, increment=0.2,
                    textvariable=self.zoom_factor_ctrl).grid(row=4, column=1, sticky="ew")

        # Iteration scaling (exponential)
        ttk.Checkbutton(sf, text="Auto-increase iterations (exponential)",
                        variable=self.auto_iter).grid(row=5, column=0, columnspan=2, sticky="w", pady=(8, 0))
        ttk.Label(sf, text="Power").grid(row=6, column=0, sticky="w", padx=(0, 6))
        ttk.Spinbox(sf, from_=0.5, to=3.0, increment=0.1,
                    textvariable=self.iter_power).grid(row=6, column=1, sticky="ew")

        # High precision
        ttk.Checkbutton(sf, text="High Precision Mode",
                        variable=self.high_prec).grid(row=7, column=0, columnspan=2, sticky="w", pady=(8, 0))
        ttk.Label(sf, text="Decimal digits").grid(row=8, column=0, sticky="w", padx=(0, 6))
        ttk.Spinbox(sf, from_=30, to=500, increment=10,
                    textvariable=self.decimal_places).grid(row=8, column=1, sticky="ew")

        # Zoom readout + actions
        ttk.Label(sf, text="Current zoom (×)").grid(row=9, column=0, sticky="w", padx=(0, 6))
        ttk.Label(sf, textvariable=self.zoom_level).grid(row=9, column=1, sticky="w")

        ttk.Button(sf, text="Generate", command=self._on_generate_click)\
            .grid(row=10, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        ttk.Button(sf, text="Reset view", command=self._reset_view)\
            .grid(row=11, column=0, columnspan=2, sticky="ew", pady=(6, 0))
        ttk.Button(sf, text="Save PNG", command=self._save_png)\
            .grid(row=12, column=0, columnspan=2, sticky="ew", pady=(6, 0))

        # Progress
        self.progress = ttk.Progressbar(self.sidebar, mode="determinate", maximum=100)
        self.progress.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(8, 0))

    # -------- Theming / drawing --------
    def _apply_dark_mode(self):
        dark = self.dark_mode.get()
        bg = "#0e1117" if dark else "#ffffff"
        fg = "#e5e7eb" if dark else "#111827"
        self.ax.set_facecolor(bg)
        self.fig.patch.set_facecolor(bg)
        for s in self.ax.spines.values():
            s.set_color(fg)
        self.ax.tick_params(colors=fg)
        self.ax.set_xticks([]); self.ax.set_yticks([])
        self.canvas.draw()

    def _plot_image(self, img):
        self.ax.clear()
        cmap = self.color_map.get() if self.enable_color.get() else "gray"
        self.ax.imshow(
            img,
            extent=[self.xmin.get(), self.xmax.get(), self.ymin.get(), self.ymax.get()],
            origin="lower",
            cmap=cmap,
            interpolation="nearest",
        )
        self.ax.set_xticks([]); self.ax.set_yticks([])
        self._apply_dark_mode()
        self.canvas.draw()
        self._ensure_selector()  # reattach selector after clear

    def _redraw_last(self):
        if self.image_data is not None:
            self._plot_image(self.image_data)

    # -------- Iteration logic --------
    def _current_iterations(self):
        """
        Exponential scaling when enabled:
            iter = base_iter * (zoom_level ** power)
        """
        zl = max(1.0, float(self.zoom_level.get()))
        base = max(10, int(self._base_iter_snapshot))
        if not self.auto_iter.get():
            return base
        power = float(self.iter_power.get())
        it = int(round(base * (zl ** power)))
        return max(10, min(it, 20000))

    # -------- Generate / Save --------
    def _on_generate_click(self):
        # Capture the user's current Max Iter as the base for auto-scaling
        self._base_iter_snapshot = int(self.max_iter.get())
        self._generate_async()

    def _generate_async(self):
        w, h = int(self.width.get()), int(self.height.get())
        xmin, xmax = float(self.xmin.get()), float(self.xmax.get())
        ymin, ymax = float(self.ymin.get()), float(self.ymax.get())
        if xmin >= xmax or ymin >= ymax:
            messagebox.showerror("Invalid bounds", "Require xmin < xmax and ymin < ymax.")
            return

        # Update zoom readout (geometric mean shrink)
        xr0, yr0 = self._start_span
        xr = xmax - xmin
        yr = ymax - ymin
        zl = ((xr0 / xr) * (yr0 / yr)) ** 0.5
        self.zoom_level.set(round(zl, 6))

        mi = self._current_iterations()

        self.status.config(text=f"Rendering {w}×{h}, {mi} iter…")
        self.progress["value"] = 0

        def work():
            t0 = time.time()
            if self.high_prec.get():
                img = mandelbrot_highprec(xmin, xmax, ymin, ymax, w, h, mi, int(self.decimal_places.get()))
            else:
                img = mandelbrot_calc(xmin, xmax, ymin, ymax, w, h, mi)
            dt = time.time() - t0
            def done():
                self.image_data = img
                self._plot_image(img)
                self.status.config(text=f"Done in {dt:.2f} s • Zoom {self.zoom_level.get()}× • Iter {mi}")
                self.progress["value"] = 100
            self.after(0, done)

        threading.Thread(target=work, daemon=True).start()

    def _save_png(self):
        if self.image_data is None:
            messagebox.showinfo("Nothing to save", "Generate an image first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG image", "*.png")])
        if not path:
            return
        buf = io.BytesIO()
        self.fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
        with open(path, "wb") as f:
            f.write(buf.getvalue())
        self.status.config(text=f"Saved PNG → {path}")

    # -------- Drag-zoom --------
    def _toggle_selector(self):
        if self.drag_zoom.get():
            self._ensure_selector()
            if self.selector:
                self.selector.set_active(True)
        else:
            if self.selector:
                self.selector.set_active(False)

    def _ensure_selector(self):
        if not self.drag_zoom.get():
            return
        # Recreate selector if axis got cleared
        if self.selector is None or self.selector.ax is not self.ax:
            self.selector = RectangleSelector(
                self.ax,
                onselect=self._on_select,
                useblit=True,
                button=[1],
                minspanx=2, minspany=2,
                spancoords='pixels',
                interactive=False,
                drag_from_anywhere=True
            )
        else:
            self.selector.set_active(True)

    def _on_select(self, eclick, erelease):
        try:
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata
            if None in (x1, y1, x2, y2):
                return
            cx = 0.5 * (x1 + x2)
            cy = 0.5 * (y1 + y2)
        except Exception:
            return

        factor = max(1.0001, float(self.zoom_factor_ctrl.get()))
        xr = self.xmax.get() - self.xmin.get()
        yr = self.ymax.get() - self.ymin.get()
        new_xr = xr / factor
        new_yr = yr / factor

        self.xmin.set(cx - new_xr / 2.0)
        self.xmax.set(cx + new_xr / 2.0)
        self.ymin.set(cy - new_yr / 2.0)
        self.ymax.set(cy + new_yr / 2.0)
        self._generate_async()

    def _reset_view(self):
        x0, x1, y0, y1 = self._defaults
        self.xmin.set(x0); self.xmax.set(x1)
        self.ymin.set(y0); self.ymax.set(y1)
        # Reset zoom level and iteration base to current Max Iter
        self.zoom_level.set(1.0)
        self._base_iter_snapshot = int(self.max_iter.get())
        self._generate_async()


if __name__ == "__main__":
    app = MandelbrotApp()
    app.mainloop()

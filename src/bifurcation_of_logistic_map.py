import time
import io
import threading
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib import transforms

# --- Constants ---
FEIGENBAUM_R = 3.569945672  # onset of chaos (accumulation point)
# A few well-known period-doubling r-values for the logistic map
KNOWN_PD = np.array([
    3.0000000000000000,  # 1 -> 2
    3.4494897430000002,  # 2 -> 4  (≈ 1 + sqrt(6))
    3.5440903590000001,  # 4 -> 8
    3.564407266,         # 8 -> 16
    3.568759400,         # ...
    3.569691610,         # ...
    FEIGENBAUM_R         # accumulation
], dtype=np.float64)


# --- Core map ---
def logistic_iter(x, r):
    return r * x * (1.0 - x)


# --- Bifurcation computation (vectorized) ---
def compute_bifurcation(r_min, r_max, num_r, burn_in, keep, x0_mode, x0_value, seed, progress_cb=None):
    rng = np.random.default_rng(seed)
    r = np.linspace(r_min, r_max, num_r, dtype=np.float64)
    x = rng.random(num_r) if x0_mode == "Random per r" else np.full(num_r, float(x0_value))
    for _ in range(burn_in):
        x = logistic_iter(x, r)
    R = np.repeat(r, keep).astype(np.float32)
    X = np.empty(num_r * keep, dtype=np.float32)
    for i in range(keep):
        x = logistic_iter(x, r)
        X[i * num_r:(i + 1) * num_r] = x.astype(np.float32)
        if progress_cb:
            progress_cb(int((i + 1) / keep * 100))
    return R, X


# --- Helpers for robust cycle detection/validation ---
def near_bifurcation(r, base=1e-6):
    """Return True if r is within an adaptive epsilon of a known bifurcation value."""
    eps = max(base, 5e-5 * (1 + abs(r)))  # tune 5e-5 for stricter/looser skip band
    return bool(np.any(np.abs(KNOWN_PD - r) < eps))


def detect_cycle_values(r, x0=0.23, burn=2000, samples=1024, tol=1e-8):
    """Return sorted approximations to the attracting cycle values for r (periodic regime).
    For r < FEIGENBAUM_R the attractor is a 1,2,4,8,... cycle. We cluster the tail
    samples to estimate the distinct levels.
    """
    # burn-in
    x = float(x0)
    for _ in range(burn):
        x = logistic_iter(x, r)
    # collect window
    vals = np.empty(samples, dtype=np.float64)
    for i in range(samples):
        x = logistic_iter(x, r)
        vals[i] = x
    vals.sort()
    if len(vals) == 0:
        return []
    # cluster close values
    clusters = []
    cur = [vals[0]]
    for v in vals[1:]:
        if abs(v - cur[-1]) <= tol:
            cur.append(v)
        else:
            clusters.append(float(np.mean(cur)))
            cur = [v]
    clusters.append(float(np.mean(cur)))
    # keep inside [0,1]
    clusters = [c for c in clusters if 0.0 <= c <= 1.0]
    # dedupe again with a slightly stricter tol
    out = []
    for c in clusters:
        if not out or abs(c - out[-1]) > tol / 5:
            out.append(c)
    return sorted(out)


def well_separated(vals, delta_sep=1e-4):
    vals = np.sort(np.asarray(vals))
    if len(vals) <= 1:
        return True
    return bool(np.all(np.diff(vals) > delta_sep))


def closes_under_iteration(r, vals, tol_close=1e-6):
    p = len(vals)
    if p == 0:
        return False
    x = float(vals[0])
    for _ in range(p):
        x = logistic_iter(x, r)
    return abs(x - vals[0]) < tol_close


def cycle_multiplier(r, pts):
    # For logistic: f'(x) = r(1-2x); multiplier over one period is product of derivatives
    pts = np.asarray(pts, dtype=np.float64)
    return float(np.prod(r * (1.0 - 2.0 * pts)))


def is_stable_cycle(r, pts, margin=1e-3):
    lam = cycle_multiplier(r, pts)
    return abs(lam) < (1.0 - margin)


# --- Tk App ---
class BifurcationApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Logistic Map — Bifurcation, Cobweb & Time Series (Tk)")
        self.geometry("1300x1000")  # default window size
        self.minsize(1200, 800)      # minimum size

        # declare Tk variables BEFORE building UI
        self.dark_mode = tk.BooleanVar(value=True)
        self.show_axes = tk.BooleanVar(value=False)
        self.show_grid = tk.BooleanVar(value=False)

        # Marker controls
        self.show_markers = tk.BooleanVar(value=False)
        self.r_mark1 = tk.DoubleVar(value=3.000)
        self.r_mark2 = tk.DoubleVar(value=3.449)

        # Time series cycle lines toggle
        self.ts_show_cycle_lines = tk.BooleanVar(value=False)

        # main layout
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        self.sidebar = ttk.Frame(self, padding=10)
        self.sidebar.grid(row=0, column=0, sticky="nsw")
        self.main = ttk.Frame(self, padding=(10, 10, 10, 10))
        self.main.grid(row=0, column=1, sticky="nsew")
        self.main.columnconfigure(0, weight=1)
        self.main.rowconfigure(0, weight=1)

        # tabs + figures
        self.tabs = ttk.Notebook(self.main)
        self.tabs.grid(row=0, column=0, sticky="nsew")
        tab1 = ttk.Frame(self.tabs); tab1.columnconfigure(0, weight=1); tab1.rowconfigure(0, weight=1)
        tab2 = ttk.Frame(self.tabs); tab2.columnconfigure(0, weight=1); tab2.rowconfigure(0, weight=1)
        tab3 = ttk.Frame(self.tabs); tab3.columnconfigure(0, weight=1); tab3.rowconfigure(0, weight=1)  # Time Series
        self.tabs.add(tab1, text="Bifurcation")
        self.tabs.add(tab2, text="Cobweb")
        self.tabs.add(tab3, text="Time Series")

        # Bifurcation fig
        self.fig_bif = Figure(figsize=(8, 6), dpi=100)
        self.ax_bif = self.fig_bif.add_subplot(111)
        self.canvas_bif = FigureCanvasTkAgg(self.fig_bif, master=tab1)
        self.canvas_bif.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        # Cobweb fig
        self.fig_cw = Figure(figsize=(8, 6), dpi=100)
        self.ax_cw = self.fig_cw.add_subplot(111)
        self.canvas_cw = FigureCanvasTkAgg(self.fig_cw, master=tab2)
        self.canvas_cw.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        # Time Series fig
        self.fig_ts = Figure(figsize=(8, 6), dpi=100)
        self.ax_ts = self.fig_ts.add_subplot(111)
        self.canvas_ts = FigureCanvasTkAgg(self.fig_ts, master=tab3)
        self.canvas_ts.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        # status
        self.status = ttk.Label(self.main, text="Ready.", anchor="w")
        self.status.grid(row=1, column=0, sticky="ew", pady=(8, 0))

        # sidebar controls (aligned grid)
        self._build_sidebar()

        # state
        self.current_R = None
        self.current_X = None

        # first render + auto-generate once
        self._apply_dark_mode()
        self._draw_bif_placeholder()
        self._draw_cobweb_placeholder()
        self._draw_timeseries_placeholder()
        self.after(50, self._generate_async)

    def _build_sidebar(self):
        # use grid: col 0 = labels, col 1 = widgets/buttons
        self.sidebar.columnconfigure(0, weight=0)
        self.sidebar.columnconfigure(1, weight=1)

        # --- Bifurcation Controls ---
        bf = ttk.LabelFrame(self.sidebar, text="Bifurcation Controls", padding=10)
        bf.grid(row=0, column=0, columnspan=2, sticky="new", pady=(0, 10))
        for c in (0, 1):
            bf.columnconfigure(c, weight=1 if c == 1 else 0)

        self.r_min = tk.DoubleVar(value=2.8)
        self.r_max = tk.DoubleVar(value=4.0)
        self.num_r = tk.IntVar(value=4000)
        self.burn_in = tk.IntVar(value=500)
        self.keep = tk.IntVar(value=1)
        self.x0_mode = tk.StringVar(value="Random per r")
        self.x0_value = tk.DoubleVar(value=0.123)
        self.seed = tk.IntVar(value=51)

        def row(r, text, widget):
            ttk.Label(bf, text=text).grid(row=r, column=0, sticky="w", pady=2, padx=(0, 6))
            widget.grid(row=r, column=1, sticky="ew", pady=2)

        row(0, "r min", ttk.Spinbox(bf, from_=0.0, to=4.0, textvariable=self.r_min, increment=0.001))
        row(1, "r max", ttk.Spinbox(bf, from_=0.0, to=4.0, textvariable=self.r_max, increment=0.001))
        row(2, "# r samples", ttk.Spinbox(bf, from_=100, to=50000, textvariable=self.num_r, increment=100))
        row(3, "Burn-in iters", ttk.Spinbox(bf, from_=0, to=10000, textvariable=self.burn_in, increment=10))
        row(4, "Points per r", ttk.Spinbox(bf, from_=1, to=1, textvariable=self.keep, increment=10))

        # initial condition (dropdown + x0 spinbox aligned)
        ic = ttk.OptionMenu(bf, self.x0_mode, self.x0_mode.get(), "Random per r", "Fixed x₀")
        row(5, "Initial condition", ic)
        row(6, "x₀ (if fixed)", ttk.Spinbox(bf, from_=0.0, to=1.0, textvariable=self.x0_value, increment=0.001))
        row(7, "Random seed", ttk.Spinbox(bf, from_=0, to=10**9, textvariable=self.seed, increment=1))

        # --- Style ---
        sf = ttk.LabelFrame(self.sidebar, text="Style", padding=10)
        sf.grid(row=1, column=0, columnspan=2, sticky="new", pady=(0, 10))
        for c in (0, 1): sf.columnconfigure(c, weight=1 if c == 1 else 0)

        self.pt_size = tk.DoubleVar(value=1.0)
        self.alpha = tk.DoubleVar(value=0.8)

        ttk.Label(sf, text="Point size").grid(row=0, column=0, sticky="w", padx=(0, 6))
        ttk.Scale(sf, from_=0.05, to=3.0, variable=self.pt_size, orient="horizontal",
                  command=lambda *_: self._redraw_bif()).grid(row=0, column=1, sticky="ew")

        ttk.Label(sf, text="Opacity").grid(row=1, column=0, sticky="w", padx=(0, 6))
        ttk.Scale(sf, from_=0.05, to=1.0, variable=self.alpha, orient="horizontal",
                  command=lambda *_: self._redraw_bif()).grid(row=1, column=1, sticky="ew")

        ttk.Checkbutton(sf, text="Dark mode", variable=self.dark_mode,
                        command=self._apply_dark_mode).grid(row=2, column=1, sticky="w", pady=(4, 0))
        ttk.Checkbutton(sf, text="Show axes", variable=self.show_axes,
                        command=self._redraw_all_plots).grid(row=3, column=1, sticky="w")
        ttk.Checkbutton(sf, text="Show grid", variable=self.show_grid,
                        command=self._redraw_all_plots).grid(row=4, column=1, sticky="w")

        # --- Markers (r-interval) ---
        mk = ttk.LabelFrame(self.sidebar, text="Markers (r-interval)", padding=10)
        mk.grid(row=2, column=0, columnspan=2, sticky="new", pady=(0, 10))
        for c in (0, 1): mk.columnconfigure(c, weight=1 if c == 1 else 0)

        ttk.Checkbutton(mk, text="Show r-interval markers", variable=self.show_markers,
                        command=self._redraw_bif).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0,4))

        ttk.Label(mk, text="r₁").grid(row=1, column=0, sticky="w", padx=(0,6))
        ttk.Spinbox(mk, from_=0.0, to=4.0, increment=0.001, textvariable=self.r_mark1,
                    command=self._redraw_bif).grid(row=1, column=1, sticky="ew", pady=2)

        ttk.Label(mk, text="r₂").grid(row=2, column=0, sticky="w", padx=(0,6))
        ttk.Spinbox(mk, from_=0.0, to=4.0, increment=0.001, textvariable=self.r_mark2,
                    command=self._redraw_bif).grid(row=2, column=1, sticky="ew", pady=2)

        # --- Cobweb ---
        cwf = ttk.LabelFrame(self.sidebar, text="Cobweb Explorer", padding=10)
        cwf.grid(row=3, column=0, columnspan=2, sticky="new", pady=(0, 10))
        for c in (0, 1): cwf.columnconfigure(c, weight=1 if c == 1 else 0)

        self.r_cw = tk.DoubleVar(value=3.7)
        self.x0_cw = tk.DoubleVar(value=0.2)
        self.steps_cw = tk.IntVar(value=60)

        def row2(r, text, widget):
            ttk.Label(cwf, text=text).grid(row=r, column=0, sticky="w", pady=2, padx=(0, 6))
            widget.grid(row=r, column=1, sticky="ew", pady=2)

        row2(0, "r (cobweb)", ttk.Spinbox(cwf, from_=0.0, to=4.0, textvariable=self.r_cw, increment=0.001,
                                           command=self._redraw_cobweb))
        row2(1, "x₀ (cobweb)", ttk.Spinbox(cwf, from_=0.0, to=1.0, textvariable=self.x0_cw, increment=0.001,
                                            command=self._redraw_cobweb))
        row2(2, "Steps", ttk.Spinbox(cwf, from_=1, to=2000, textvariable=self.steps_cw, increment=1,
                                     command=self._redraw_cobweb))

        # --- Time Series ---
        tsf = ttk.LabelFrame(self.sidebar, text="Time Series (xₙ vs n)", padding=10)
        tsf.grid(row=4, column=0, columnspan=2, sticky="new", pady=(0, 10))
        for c in (0, 1): tsf.columnconfigure(c, weight=1 if c == 1 else 0)

        self.r_ts = tk.DoubleVar(value=3.7)
        self.x0_ts = tk.DoubleVar(value=0.2)
        self.steps_ts = tk.IntVar(value=100)

        ttk.Label(tsf, text="r (time series)").grid(row=0, column=0, sticky="w", padx=(0,6))
        ttk.Spinbox(tsf, from_=0.0, to=4.0, textvariable=self.r_ts, increment=0.001,
                    command=self._redraw_timeseries).grid(row=0, column=1, sticky="ew", pady=2)
        ttk.Label(tsf, text="x₀ (time series)").grid(row=1, column=0, sticky="w", padx=(0,6))
        ttk.Spinbox(tsf, from_=0.0, to=1.0, textvariable=self.x0_ts, increment=0.001,
                    command=self._redraw_timeseries).grid(row=1, column=1, sticky="ew", pady=2)
        ttk.Label(tsf, text="Steps n").grid(row=2, column=0, sticky="w", padx=(0,6))
        ttk.Spinbox(tsf, from_=1, to=10000, textvariable=self.steps_ts, increment=1,
                    command=self._redraw_timeseries).grid(row=2, column=1, sticky="ew", pady=2)

        ttk.Checkbutton(tsf, text="Show cycle levels (red) for periodic r",
                        variable=self.ts_show_cycle_lines,
                        command=self._redraw_timeseries).grid(row=3, column=0, columnspan=2, sticky="w", pady=(6,0))

        # --- Actions (aligned with inputs) ---
        ttk.Button(self.sidebar, text="Generate", command=self._generate_async)\
            .grid(row=5, column=1, sticky="ew")
        ttk.Button(self.sidebar, text="Save PNG (bifurcation)", command=self._save_png_bif)\
            .grid(row=6, column=1, sticky="ew", pady=(6, 0))
        ttk.Button(self.sidebar, text="Save PNG (cobweb)", command=self._save_png_cw)\
            .grid(row=7, column=1, sticky="ew", pady=(6, 0))
        ttk.Button(self.sidebar, text="Save PNG (time series)", command=self._save_png_ts)\
            .grid(row=8, column=1, sticky="ew", pady=(6, 0))

        # progress + tips
        self.progress = ttk.Progressbar(self.sidebar, mode="determinate", maximum=100)
        self.progress.grid(row=9, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        ttk.Label(self.sidebar,
                  text="Tip: total points = samples × points/r.\nKeep under ~1M for smooth plotting.\nRandom x₀ reveals full structure.",
                  justify="left").grid(row=10, column=0, columnspan=2, sticky="ew", pady=(8, 0))

    # --- Theming / drawing helpers ---
    def _apply_dark_mode(self):
        dark = self.dark_mode.get()
        bg = "#0e1117" if dark else "#ffffff"
        fg = "#e5e7eb" if dark else "#111827"
        for ax, fig in [(self.ax_bif, self.fig_bif), (self.ax_cw, self.fig_cw), (self.ax_ts, self.fig_ts)]:
            ax.set_facecolor(bg); fig.patch.set_facecolor(bg)
            for s in ax.spines.values(): s.set_color(fg)
            ax.tick_params(colors=fg)
            ax.yaxis.label.set_color(fg); ax.xaxis.label.set_color(fg); ax.title.set_color(fg)
        # force redraw
        self.canvas_bif.draw()
        self.canvas_cw.draw()
        self.canvas_ts.draw()

    def _toggle_axes_grid(self, ax):
        if self.show_axes.get():
            for s in ax.spines.values(): s.set_visible(True)
            ax.set_xticks(ax.get_xticks()); ax.set_yticks(ax.get_yticks())
        else:
            for s in ax.spines.values(): s.set_visible(False)
            ax.set_xticks([]); ax.set_yticks([])
        ax.grid(self.show_grid.get(), linewidth=0.5, alpha=0.3)

    def _draw_bif_placeholder(self):
        self.ax_bif.clear()
        self.ax_bif.set_title("Bifurcation Diagram")
        self.ax_bif.set_xlabel("r"); self.ax_bif.set_ylabel("xₙ")
        self.ax_bif.set_xlim(2.8, 4.0); self.ax_bif.set_ylim(0, 1)
        self._toggle_axes_grid(self.ax_bif)
        self._apply_dark_mode()
        self.canvas_bif.draw()

    def _draw_cobweb_placeholder(self):
        self.ax_cw.clear()
        self.ax_cw.set_title("Cobweb Plot")
        self.ax_cw.set_xlabel("x"); self.ax_cw.set_ylabel("f(x)")
        self.ax_cw.set_xlim(0, 1); self.ax_cw.set_ylim(0, 1)
        self._toggle_axes_grid(self.ax_cw)
        self._apply_dark_mode()
        self.canvas_cw.draw()

    def _draw_timeseries_placeholder(self):
        self.ax_ts.clear()
        self.ax_ts.set_title("Time Series")
        self.ax_ts.set_xlabel("n"); self.ax_ts.set_ylabel("xₙ")
        self.ax_ts.set_xlim(0, max(10, int(self.steps_ts.get())))
        self.ax_ts.set_ylim(0, 1)
        self._toggle_axes_grid(self.ax_ts)
        self._apply_dark_mode()
        self.canvas_ts.draw()

    def _redraw_bif(self):
        if self.current_R is None:
            return
        self._plot_bif(self.current_R, self.current_X)

    def _plot_bif(self, R, X):
        self.ax_bif.clear()
        # high-contrast point color
        point_color = "white" if self.dark_mode.get() else "black"
        self.ax_bif.scatter(R, X, s=self.pt_size.get(), alpha=self.alpha.get(),
                            c=point_color, marker=".", linewidths=0)
        self.ax_bif.set_xlim(float(self.r_min.get()), float(self.r_max.get()))
        self.ax_bif.set_ylim(0, 1)
        self.ax_bif.set_title("Bifurcation Diagram of the Logistic Map")
        self.ax_bif.set_xlabel("r")
        self.ax_bif.set_ylabel("xₙ")

        # --- r-interval markers + labels below axis ---
        if self.show_markers.get():
            r1 = float(self.r_mark1.get())
            r2 = float(self.r_mark2.get())

            # draw vertical lines
            self.ax_bif.axvline(r1, color="red", linewidth=1.2, alpha=0.9, linestyle=":")
            self.ax_bif.axvline(r2, color="red", linewidth=1.2, alpha=0.9, linestyle=":")

            # label color matches tick color for theme
            label_color = "#e5e7eb" if self.dark_mode.get() else "#111827"

            # blended transform: x in data coords, y in axes coords
            trans = transforms.blended_transform_factory(self.ax_bif.transData, self.ax_bif.transAxes)
            # place labels slightly below the axis (y = -0.08)
            for rv in (r1, r2):
                self.ax_bif.text(rv, -0.08, f"{rv:.3f}", transform=trans,
                                 ha="center", va="top", color=label_color, fontsize=9, clip_on=False)

        self._toggle_axes_grid(self.ax_bif)
        self._apply_dark_mode()
        self.canvas_bif.draw()  # <- force draw

    def _redraw_cobweb(self, *_):
        r = float(self.r_cw.get()); x0 = float(self.x0_cw.get()); steps = int(self.steps_cw.get())
        self.ax_cw.clear()
        xs = np.linspace(0, 1, 1000)
        fx = logistic_iter(xs, r)
        line_color = "white" if self.dark_mode.get() else "black"
        self.ax_cw.plot(xs, fx, linewidth=1.5, color=line_color)
        self.ax_cw.plot(xs, xs, linewidth=1.0, color=line_color)

        x = x0
        for _ in range(steps):
            y = logistic_iter(x, r)
            self.ax_cw.plot([x, x], [x, y], linewidth=0.8, color=line_color)
            self.ax_cw.plot([x, y], [y, y], linewidth=0.8, color=line_color)
            x = y

        self.ax_cw.set_xlim(0, 1); self.ax_cw.set_ylim(0, 1)
        self.ax_cw.set_xlabel("x"); self.ax_cw.set_ylabel("f(x)")
        self.ax_cw.set_title(f"Cobweb Plot (r = {r:.4f}, x₀ = {x0:.3f}, steps = {steps})")
        self._toggle_axes_grid(self.ax_cw)
        self._apply_dark_mode()
        self.canvas_cw.draw()   # <- force draw

    def _redraw_timeseries(self, *_):
        r = float(self.r_ts.get()); x0 = float(self.x0_ts.get()); steps = int(self.steps_ts.get())
        # compute trajectory x_n
        xs = np.empty(steps + 1, dtype=np.float64)
        xs[0] = x0
        for n in range(steps):
            xs[n + 1] = logistic_iter(xs[n], r)
        n_vals = np.arange(steps + 1)

        self.ax_ts.clear()
        line_color = "white" if self.dark_mode.get() else "black"
        self.ax_ts.plot(n_vals, xs, linewidth=1.5, color=line_color)
        self.ax_ts.scatter(n_vals, xs, s=10, color=line_color, alpha=0.7)
        self.ax_ts.set_xlim(0, max(10, steps))
        self.ax_ts.set_ylim(0, 1)
        self.ax_ts.set_xlabel("n")
        self.ax_ts.set_ylabel("xₙ")
        self.ax_ts.set_title(f"Time Series (r = {r:.4f}, x₀ = {x0:.3f}, steps = {steps})")

        # draw cycle levels if requested and robustly away from bifurcation/chaos
        if self.ts_show_cycle_lines.get():
            if r < FEIGENBAUM_R and not near_bifurcation(r):
                cycle_vals = detect_cycle_values(r, x0=x0)
                if (cycle_vals and well_separated(cycle_vals, 1e-4)
                        and closes_under_iteration(r, cycle_vals, 1e-6)
                        and is_stable_cycle(r, cycle_vals, margin=1e-3)):
                    for y in cycle_vals:
                        # horizontal dashed line
                        self.ax_ts.axhline(y, xmin=0, xmax=1, linestyle='--', linewidth=1.2, color='red', alpha=0.9)
                        # label slightly left of axis: x in axes coords, y in data coords
                        trans = transforms.blended_transform_factory(self.ax_ts.transAxes, self.ax_ts.transData)
                        self.ax_ts.text(-0.08, y, f"{y:.4f}", transform=trans, ha='right', va='center',
                                        color='red', fontsize=9, clip_on=False)
                    self.status.config(text=f"Cycle levels: {len(cycle_vals)} value(s) shown.")
                else:
                    self.status.config(text="Cycle not drawn (near boundary / not clean / not stable).")
            else:
                self.status.config(text="Skipped cycle lines near bifurcation or in chaos.")

        self._toggle_axes_grid(self.ax_ts)
        self._apply_dark_mode()
        self.canvas_ts.draw()

    def _redraw_all_plots(self):
        self._redraw_bif()
        self._redraw_cobweb()
        self._redraw_timeseries()

    # --- actions ---
    def _generate_async(self):
        try:
            rmin = float(self.r_min.get()); rmax = float(self.r_max.get())
            assert 0 <= rmin < rmax <= 4.0
        except Exception:
            messagebox.showerror("Invalid input", "Check r min/max (0 ≤ r_min < r_max ≤ 4).")
            return
        num_r = int(self.num_r.get()); burn_in = int(self.burn_in.get()); keep = int(self.keep.get())
        x0_mode = self.x0_mode.get(); x0_value = float(self.x0_value.get()); seed = int(self.seed.get())
        total_pts = num_r * keep
        self.status.config(text=f"Computing ~{total_pts:,} points…"); self.progress["value"] = 0

        def progress_cb(p):
            self.after(0, lambda: self.progress.configure(value=p))

        def work():
            t0 = time.time()
            try:
                R, X = compute_bifurcation(rmin, rmax, num_r, burn_in, keep, x0_mode, x0_value, seed, progress_cb)
            except Exception as e:
                self.after(0, lambda: messagebox.showerror("Error", str(e)))
                return
            dt = time.time() - t0

            def done():
                self.current_R, self.current_X = R, X
                self._plot_bif(R, X)
                self._redraw_cobweb()
                self._redraw_timeseries()
                self.status.config(text=f"Done in {dt:.2f} s • Points: {total_pts:,}")
                self.progress["value"] = 100

            self.after(0, done)

        threading.Thread(target=work, daemon=True).start()

    def _save_png_bif(self):
        if self.current_R is None:
            messagebox.showinfo("Nothing to save", "Generate the bifurcation diagram first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG image", "*.png")])
        if not path:
            return
        buf = io.BytesIO(); self.fig_bif.savefig(buf, format="png", dpi=200, bbox_inches="tight")
        with open(path, "wb") as f:
            f.write(buf.getvalue())
        self.status.config(text=f"Saved bifurcation PNG → {path}")

    def _save_png_cw(self):
        path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG image", "*.png")])
        if not path:
            return
        buf = io.BytesIO(); self.fig_cw.savefig(buf, format="png", dpi=200, bbox_inches="tight")
        with open(path, "wb") as f:
            f.write(buf.getvalue())
        self.status.config(text=f"Saved cobweb PNG → {path}")

    def _save_png_ts(self):
        path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG image", "*.png")])
        if not path:
            return
        buf = io.BytesIO(); self.fig_ts.savefig(buf, format="png", dpi=200, bbox_inches="tight")
        with open(path, "wb") as f:
            f.write(buf.getvalue())
        self.status.config(text=f"Saved time series PNG → {path}")


if __name__ == "__main__":
    app = BifurcationApp()
    app.mainloop()

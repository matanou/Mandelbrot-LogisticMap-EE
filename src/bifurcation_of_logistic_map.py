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


def logistic_iter(x, r):
    return r * x * (1.0 - x)


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


class BifurcationApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Logistic Map — Bifurcation & Cobweb (Tk)")
        self.geometry("1200x800")
        self.minsize(1000, 700)

        # declare Tk variables BEFORE building UI
        self.dark_mode = tk.BooleanVar(value=True)
        self.show_axes = tk.BooleanVar(value=False)
        self.show_grid = tk.BooleanVar(value=False)

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
        self.tabs.add(tab1, text="Bifurcation")
        self.tabs.add(tab2, text="Cobweb")

        self.fig_bif = Figure(figsize=(7, 5), dpi=100)
        self.ax_bif = self.fig_bif.add_subplot(111)
        self.canvas_bif = FigureCanvasTkAgg(self.fig_bif, master=tab1)
        self.canvas_bif.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        self.fig_cw = Figure(figsize=(6, 5), dpi=100)
        self.ax_cw = self.fig_cw.add_subplot(111)
        self.canvas_cw = FigureCanvasTkAgg(self.fig_cw, master=tab2)
        self.canvas_cw.get_tk_widget().grid(row=0, column=0, sticky="nsew")

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
        row(3, "Burn-in iters", ttk.Spinbox(bf, from_=0, to=10000, textvariable=self.burn_in, increment=50))
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

        self.pt_size = tk.DoubleVar(value=1.0)   # more visible defaults
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
                        command=self._redraw_bif).grid(row=3, column=1, sticky="w")
        ttk.Checkbutton(sf, text="Show grid", variable=self.show_grid,
                        command=self._redraw_bif).grid(row=4, column=1, sticky="w")

        # --- Cobweb ---
        cwf = ttk.LabelFrame(self.sidebar, text="Cobweb Explorer", padding=10)
        cwf.grid(row=2, column=0, columnspan=2, sticky="new", pady=(0, 10))
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

        # --- Actions (aligned with inputs) ---
        ttk.Button(self.sidebar, text="Generate", command=self._generate_async)\
            .grid(row=3, column=1, sticky="ew")
        ttk.Button(self.sidebar, text="Save PNG (bifurcation)", command=self._save_png_bif)\
            .grid(row=4, column=1, sticky="ew", pady=(6, 0))
        ttk.Button(self.sidebar, text="Save PNG (cobweb)", command=self._save_png_cw)\
            .grid(row=5, column=1, sticky="ew", pady=(6, 0))

        # progress + tips
        self.progress = ttk.Progressbar(self.sidebar, mode="determinate", maximum=100)
        self.progress.grid(row=6, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        ttk.Label(self.sidebar,
                  text="Tip: total points = samples × points/r.\nKeep under ~1M for smooth plotting.\nRandom x₀ reveals full structure.",
                  justify="left").grid(row=7, column=0, columnspan=2, sticky="ew", pady=(8, 0))

    # --- Theming / drawing helpers ---
    def _apply_dark_mode(self):
        dark = self.dark_mode.get()
        bg = "#0e1117" if dark else "#ffffff"
        fg = "#e5e7eb" if dark else "#111827"
        for ax, fig in [(self.ax_bif, self.fig_bif), (self.ax_cw, self.fig_cw)]:
            ax.set_facecolor(bg); fig.patch.set_facecolor(bg)
            for s in ax.spines.values(): s.set_color(fg)
            ax.tick_params(colors=fg)
            ax.yaxis.label.set_color(fg); ax.xaxis.label.set_color(fg); ax.title.set_color(fg)
        # force redraw
        self.canvas_bif.draw()
        self.canvas_cw.draw()

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

    def _redraw_bif(self):
        if self.current_R is None: return
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
        self.ax_bif.set_xlabel("r"); self.ax_bif.set_ylabel("xₙ")
        self._toggle_axes_grid(self.ax_bif)
        self._apply_dark_mode()
        self.canvas_bif.draw()   # <- force draw

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

    # --- actions ---
    def _generate_async(self):
        try:
            rmin = float(self.r_min.get()); rmax = float(self.r_max.get())
            assert 0 <= rmin < rmax <= 4.0
        except Exception:
            messagebox.showerror("Invalid input", "Check r min/max (0 ≤ r_min < r_max ≤ 4)."); return
        num_r = int(self.num_r.get()); burn_in = int(self.burn_in.get()); keep = int(self.keep.get())
        x0_mode = self.x0_mode.get(); x0_value = float(self.x0_value.get()); seed = int(self.seed.get())
        total_pts = num_r * keep
        self.status.config(text=f"Computing ~{total_pts:,} points…"); self.progress["value"] = 0

        def progress_cb(p): self.after(0, lambda: self.progress.configure(value=p))

        def work():
            t0 = time.time()
            try:
                R, X = compute_bifurcation(rmin, rmax, num_r, burn_in, keep, x0_mode, x0_value, seed, progress_cb)
            except Exception as e:
                self.after(0, lambda: messagebox.showerror("Error", str(e))); return
            dt = time.time() - t0
            def done():
                self.current_R, self.current_X = R, X
                self._plot_bif(R, X)
                self.status.config(text=f"Done in {dt:.2f} s • Points: {total_pts:,}")
                self.progress["value"] = 100
            self.after(0, done)

        threading.Thread(target=work, daemon=True).start()

    def _save_png_bif(self):
        if self.current_R is None:
            messagebox.showinfo("Nothing to save", "Generate the bifurcation diagram first."); return
        path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG image", "*.png")])
        if not path: return
        buf = io.BytesIO(); self.fig_bif.savefig(buf, format="png", dpi=200, bbox_inches="tight")
        with open(path, "wb") as f: f.write(buf.getvalue())
        self.status.config(text=f"Saved bifurcation PNG → {path}")

    def _save_png_cw(self):
        path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG image", "*.png")])
        if not path: return
        buf = io.BytesIO(); self.fig_cw.savefig(buf, format="png", dpi=200, bbox_inches="tight")
        with open(path, "wb") as f: f.write(buf.getvalue())
        self.status.config(text=f"Saved cobweb PNG → {path}")


if __name__ == "__main__":
    app = BifurcationApp()
    app.mainloop()

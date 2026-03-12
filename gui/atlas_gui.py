#!/usr/bin/env python3
"""
ATLAS-GP  —  Universal Cluster Structure Optimizer
Session 5: Production GUI

Layout: left sidebar (all input) + right tabbed area (Monitor | Results | Log)
Follows scientific structure-search software conventions (USPEX/CALYPSO/ASE style).

Every CLI flag has exactly one corresponding widget.
JSON protocol wired precisely to main.cpp output:
  header  : {type, formula, N, potential, runs, gens}
  progress: {type, run, gen, total, best, basins}
  result  : {type, formula, potential, best_energy, n_isomers, elapsed}
XYZ comment: "Isomer=N Energy=E Formula=F"
"""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import subprocess, threading, json, os, sys, shutil, argparse, time
from pathlib import Path

CONFIG_FILE = Path.home() / ".atlas_gp_config.json"

# ── Palette: dark scientific instrument aesthetic ─────────────────────────────
C = {
    "bg":   "#0f1117",   # deepest background
    "bg1":  "#171c26",   # sidebar
    "bg2":  "#1e2433",   # panels / cards
    "bg3":  "#2a3040",   # input fields
    "bg4":  "#363d50",   # borders / dividers
    "fg":   "#d4dbe8",   # primary text
    "fg2":  "#7d8ba3",   # labels
    "fg3":  "#45506a",   # hints / disabled
    "blue": "#4d9fff",   # accent
    "blue2":"#2563eb",   # button fill
    "grn":  "#36c97a",   # success / global min
    "yel":  "#e8b84b",   # warning / running
    "red":  "#e05555",   # error / stop
    "pur":  "#a78bfa",   # basins
    "ora":  "#f59e40",   # secondary
    "tel":  "#2dd4bf",   # tertiary
}
TRACES = ["#4d9fff","#36c97a","#f59e40","#a78bfa","#2dd4bf","#e05555",
          "#fbbf24","#fb7185","#34d399","#60a5fa","#c084fc","#f97316"]

FUI  = ("Segoe UI", 9)
FBLD = ("Segoe UI", 9, "bold")
FH1  = ("Segoe UI", 11, "bold")
FMONO= ("Consolas", 9)
FBIG = ("Consolas", 13, "bold")

ELEM_COL = {"H":"#ffffff","C":"#909090","N":"#3050f8","O":"#ff0d0d",
            "F":"#90e050","P":"#ff8000","S":"#ffff30","Cl":"#1ff01f",
            "Fe":"#e06633","Ni":"#50d050","Cu":"#c88033","Ag":"#c0c0c0",
            "Au":"#ffd123","Pt":"#d0d0e0","Ar":"#80d1e3","Si":"#f0c8a0"}
def ecol(sym): return ELEM_COL.get(sym, C["blue"])

def mk_entry(parent, var, w=9, big=False):
    fn = FBIG if big else FMONO
    return tk.Entry(parent, textvariable=var, bg=C["bg3"], fg=C["fg"],
                    insertbackground=C["fg"], relief="flat", bd=0,
                    highlightthickness=1, highlightbackground=C["bg4"],
                    highlightcolor=C["blue"], font=fn, width=w)

def mk_lbl(parent, text, fg=None, font=None, bg=None, **kw):
    return tk.Label(parent, text=text, bg=bg or C["bg1"],
                    fg=fg or C["fg2"], font=font or FUI, **kw)

def mk_chk(parent, text, var, bg=None):
    b = bg or C["bg1"]
    return tk.Checkbutton(parent, text=text, variable=var, bg=b,
                          fg=C["fg2"], selectcolor=C["bg3"],
                          activebackground=b, activeforeground=C["fg"],
                          font=FUI, relief="flat", bd=0, highlightthickness=0)

def mk_btn(parent, text, cmd, bg=None, fg=None, **kw):
    font = kw.pop("font", FUI)   # allow callers to override font
    return tk.Button(parent, text=text, command=cmd,
                     bg=bg or C["bg3"], fg=fg or C["fg"], font=font,
                     relief="flat", bd=0, cursor="hand2",
                     activebackground=C["bg4"], activeforeground=C["fg"],
                     highlightthickness=0, **kw)

def divider(parent, text, bg=None):
    bg = bg or C["bg1"]
    f = tk.Frame(parent, bg=bg)
    tk.Label(f, text=text, fg=C["blue"], bg=bg, font=FBLD).pack(side="left", padx=(0,6))
    tk.Frame(f, bg=C["bg4"], height=1).pack(side="left", fill="x", expand=True, pady=5)
    return f

def param_row(parent, label, var, w=9, tip="", bg=None, registry=None):
    bg = bg or C["bg1"]
    r = tk.Frame(parent, bg=bg)
    tk.Label(r, text=label, fg=C["fg2"], bg=bg, font=FUI,
             width=20, anchor="w").pack(side="left", padx=(0,4))
    e = mk_entry(r, var, w=w)
    e.pack(side="left", padx=(0,6))
    if registry is not None:
        registry[id(var)] = e     # map id(StringVar) → Entry widget
    if tip:
        tk.Label(r, text=tip, fg=C["fg3"], bg=bg,
                 font=("Segoe UI", 8)).pack(side="left")
    return r


# ── Tooltip (hover-over help) ────────────────────────────────────────────────
class Tooltip:
    """Lightweight tooltip shown on widget hover."""
    def __init__(self, widget, text):
        self.widget = widget
        self.text   = text
        self._tip   = None
        widget.bind("<Enter>", self._show)
        widget.bind("<Leave>", self._hide)
        widget.bind("<Destroy>", self._hide)

    def _show(self, _=None):
        if self._tip:
            return
        x = self.widget.winfo_rootx() + 20
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        self._tip.wm_geometry(f"+{x}+{y}")
        self._tip.configure(bg=C["bg4"])
        tk.Label(self._tip, text=self.text, bg=C["bg2"], fg=C["fg"],
                 font=("Segoe UI", 8), relief="flat", bd=0,
                 padx=8, pady=4, justify="left",
                 wraplength=320).pack()

    def _hide(self, _=None):
        if self._tip:
            self._tip.destroy()
            self._tip = None


def tip(widget, text):
    """Attach a tooltip to widget; returns widget for chaining."""
    Tooltip(widget, text)
    return widget


# ── Parameter help strings ───────────────────────────────────────────────────
HELP = {
    "formula":       "Chemical formula of the cluster to optimise.\n"
                     "Examples: Au13  Ni55  Au32Ag6  Fe12O18  Pt13\n"
                     "Any element from H to Og is supported.",
    "potential":     "Energy model used to evaluate structures.\n"
                     "AUTO selects the best available for your formula:\n"
                     "  LJ    → noble gases (Ar, Ne, Kr)\n"
                     "  MORSE → diatomics, simple metals\n"
                     "  GUPTA → transition metals (Au, Ni, Cu, Pt…)\n"
                     "  UFF   → mixed-element / organic clusters\n"
                     "  MACE  → universal ML potential (most accurate, slower)",
    "runs":          "Number of independent basin-hopping runs.\n"
                     "Each starts from a different random structure.\n"
                     "More runs = better coverage, but longer wall-time.\n"
                     "Recommended: 8–32 for small clusters (N < 20),\n"
                     "4–8 for large clusters (N > 50).",
    "gens":          "Generations (basin-hopping steps) per run.\n"
                     "Each generation: propose → minimise → accept/reject.\n"
                     "300 is a good default; use 1000 for hard landscapes\n"
                     "(e.g. LJ38, LJ75, Au55 amorphous↔crystalline).",
    "workers":       "Number of parallel CPU threads (OpenMP).\n"
                     "Set to the number of physical cores on your machine.\n"
                     "Default 4 is conservative; 8–16 is common.",
    "seed":          "Random number seed for reproducibility.\n"
                     "0 = random seed from system clock.\n"
                     "Set a fixed integer to reproduce a previous run exactly.",
    "bht":           "Basin-hopping Metropolis acceptance temperature.\n"
                     "Higher → accept uphill moves more readily (more exploration).\n"
                     "Lower  → greedy descent (faster but gets trapped).\n"
                     "Range: 0.1 – 5.0.  Default 1.0 in reduced units.",
    "bhd":           "Multiplicative cooling factor applied each generation.\n"
                     "T(g) = T(0) × decay^g.  0.9998 cools very slowly.\n"
                     "Decrease (e.g. 0.995) for shorter runs.",
    "step":          "Perturbation step size as a fraction of mean σ.\n"
                     "Larger → bigger moves, less refined minima.\n"
                     "Range: 0.05 – 0.80.  Default 0.25.",
    "nfrac":         "Fraction of atoms displaced per perturbation step.\n"
                     "0.35 displaces ~35% of atoms each step.\n"
                     "Lower = more local refinement; higher = more diffusion.",
    "grprob":        "Probability of using bond-topology (graph) moves\n"
                     "instead of coordinate perturbation.\n"
                     "Graph moves rearrange bonded fragments; useful for\n"
                     "molecular or heteroatomic clusters.",
    "bondcut":       "Bond cutoff for graph connectivity, as a multiple of σ.\n"
                     "Atoms within bondcut × σ are considered bonded.\n"
                     "Default 1.35 works for most systems.",
    "tablam":        "Taboo memory repulsion strength λ.\n"
                     "Adds λ × exp(–d²/2δ²) to effective energy of\n"
                     "structures similar to recently visited ones.\n"
                     "Higher λ → stronger push away from visited regions.",
    "tabdel":        "Taboo memory width δ (in fingerprint distance units).\n"
                     "Controls how broadly the repulsion spreads.\n"
                     "Larger δ → repulsion affects more similar structures.",
    "kappa":         "GP LCB acquisition exploration/exploit balance κ.\n"
                     "LCB = μ – κ·σ  (lower = better in minimisation).\n"
                     "Higher κ → more exploration of uncertain regions.\n"
                     "Range: 0.5 – 5.0.  Default 2.0.",
    "cands":         "Number of random candidate moves scored by GP per step.\n"
                     "Higher → better GP-guided proposals, but slower.\n"
                     "Range: 2 – 20.  Default 5.",
    "popsize":       "Size of the shared elite population for crossover.\n"
                     "Crossover takes fragments from two population members\n"
                     "and recombines them (Session 4 innovation).\n"
                     "Range: 5 – 50.  Default 20.",
    "rmsd":          "RMSD threshold for isomer uniqueness (×mean σ).\n"
                     "Structures with Kabsch RMSD < threshold × σ\n"
                     "are considered identical and deduplicated.\n"
                     "Lower = finer discrimination; higher = fewer isomers.",
    "maxiso":        "Maximum number of distinct isomers to report.\n"
                     "The top-N by energy within the energy window are kept.\n"
                     "Default 30.",
    "ewin":          "Energy window above global minimum.\n"
                     "Isomers with E > E_min + window are discarded.\n"
                     "Expressed as a fraction (0.40 = 40%).",
    "ard":           "Automatic Relevance Determination (Session 4).\n"
                     "Learns which descriptor dimensions are most informative\n"
                     "for predicting energy, and weights the GP kernel accordingly.\n"
                     "Strongly recommended: leave enabled.",
    "ngd":           "Natural Gradient Descent for GP hyperparameters (Session 4).\n"
                     "Replaces grid search with gradient-based optimisation\n"
                     "in Riemannian (distribution) space — 2–3× faster convergence.\n"
                     "Leave enabled unless debugging.",
    "pop":           "Population crossover — evolutionary layer (Session 4).\n"
                     "Maintains a pool of elite structures; periodically combines\n"
                     "fragments from two parents to create children.\n"
                     "Proven to find LJ38/LJ75 global minima much faster.",
    "homotop":       "Homotop search for heteroatomic clusters.\n"
                     "Swaps chemical identities of atoms (not positions)\n"
                     "to find the optimal chemical ordering.\n"
                     "Only active when formula has > 1 element.",
    "eigenvec":      "Lanczos eigenvector-following escape (Session 3).\n"
                     "Computes the softest vibrational mode of the current\n"
                     "local minimum and escapes along it.\n"
                     "~40% lower barriers crossed vs random kicks.",
    "pt":            "Parallel Tempering — multiple replicas at different\n"
                     "temperatures exchange configurations (Session 3).\n"
                     "Hot replicas explore; cold replicas refine.\n"
                     "Provably samples Boltzmann distribution correctly.",
    "ptreplicas":    "Number of temperature replicas per PT group.\n"
                     "Range: 2 – 8.  Default 4.\n"
                     "More replicas = better mixing, more compute.",
    "ptmin":         "Coldest PT replica temperature (reduced units).\n"
                     "Should be ≤ BH temperature.  Default 0.3.",
    "ptmax":         "Hottest PT replica temperature (reduced units).\n"
                     "Should be high enough to escape barriers.  Default 5.0.",
    "active":        "Active learning — MACE oracle refinement (Session 3).\n"
                     "When GP uncertainty is high, evaluates structure with\n"
                     "MACE (accurate ML potential) instead of fast potential.\n"
                     "Requires MACE to be installed.",
    "oracle":        "Maximum number of MACE oracle calls.\n"
                     "Each call is expensive (~1–10 s depending on cluster size).\n"
                     "Budget is allocated to highest-uncertainty structures.",
    "outdir":        "Directory where results are written:\n"
                     "  <formula>_isomers.xyz   — all distinct isomers\n"
                     "  <formula>_global_min.xyz — best structure\n"
                     "Created automatically if it does not exist.",
}


class ATLASApp(tk.Tk):

    def __init__(self, binary, mace_script):
        super().__init__()
        self.binary = binary
        self.mace_script = mace_script
        self.mace_python = sys.executable  # updated if py -3.11 is used
        self.proc = None
        self.running = False
        self.isomers = []
        self.best_per_run = {}
        self.run_total = {}
        self.global_min = None
        self.n_runs = 8
        self.total_gens = 300
        # Protocol additions
        self._run_start_time = None   # for ETA
        self._run_history    = []     # list of {formula, best_e, n_iso, elapsed, pot}
        self._entry_reg      = {}     # StringVar → Entry widget, for inline error highlight
        self._define_vars()
        self._load_config()          # restore last session's settings
        self.title("ATLAS-GP  —  Universal Cluster Optimizer")
        self.configure(bg=C["bg"])
        self.geometry("1200x800")
        self.minsize(960, 640)
        self._build_ui()
        self.after(400, self._startup)

    def _define_vars(self):
        S, B = tk.StringVar, tk.BooleanVar
        # core
        self.vf   = S(value="Au13");       self.vp   = S(value="AUTO")
        # search
        self.vruns= S(value="8");          self.vgens= S(value="300")
        self.vwrk = S(value="4");          self.vseed= S(value="0")
        # basin-hopping
        self.vbht = S(value="1.0");        self.vbhd = S(value="0.9998")
        # perturbation
        self.vstp = S(value="0.25");       self.vnfr = S(value="0.35")
        self.vgrp = S(value="0.30");       self.vbnd = S(value="1.35")
        # taboo
        self.vtlm = S(value="0.05");       self.vtld = S(value="0.25")
        # gp
        self.bgp  = B(value=True)
        self.vkap = S(value="2.0");        self.vcnd = S(value="5")
        self.bard = B(value=True);         self.bngd = B(value=True)
        # population
        self.bpop = B(value=True);         self.vpop = S(value="20")
        # filter
        self.vrms = S(value="0.20");       self.vmxi = S(value="30")
        self.vewn = S(value="0.40")
        # features
        self.bhtop= B(value=True);         self.beig = B(value=True)
        # parallel tempering
        self.bpt  = B(value=False)
        self.vptrp= S(value="4");          self.vptlo= S(value="0.3")
        self.vpthi= S(value="5.0")
        # active learning
        self.bact = B(value=False);        self.vorc = S(value="50")
        # output
        self.vout = S(value="atlas_results"); self.bvrb= B(value=False)

    def _build_ui(self):
        # Menubar
        mb = tk.Menu(self, bg=C["bg1"], fg=C["fg"], activebackground=C["blue2"],
                     activeforeground="white", relief="flat", bd=0)
        self.config(menu=mb)

        fm = tk.Menu(mb, tearoff=0, bg=C["bg1"], fg=C["fg"],
                     activebackground=C["blue2"], activeforeground="white")
        mb.add_cascade(label="File", menu=fm)
        fm.add_command(label="Save config…    Ctrl+S", command=self._save_config_as)
        fm.add_command(label="Load config…    Ctrl+O", command=self._load_config_from)
        fm.add_separator()
        fm.add_command(label="Run history     Ctrl+H", command=self._history_window)
        fm.add_separator()
        fm.add_command(label="Exit", command=self._on_close)

        pm = tk.Menu(mb, tearoff=0, bg=C["bg1"], fg=C["fg"],
                     activebackground=C["blue2"], activeforeground="white")
        mb.add_cascade(label="Profiles", menu=pm)
        for name in self.PRESETS:
            pm.add_command(label=name, command=lambda n=name: self._apply_preset(n))

        hm = tk.Menu(mb, tearoff=0, bg=C["bg1"], fg=C["fg"],
                     activebackground=C["blue2"], activeforeground="white")
        mb.add_cascade(label="Help", menu=hm)
        hm.add_command(label="Keyboard shortcuts", command=self._show_shortcuts)

        self._bind_shortcuts()

        # Header
        hdr = tk.Frame(self, bg=C["bg1"], height=44)
        hdr.pack(fill="x"); hdr.pack_propagate(False)
        tk.Label(hdr, text="ATLAS-GP", fg=C["blue"], bg=C["bg1"],
                 font=("Consolas",14,"bold")).pack(side="left", padx=(14,0), pady=8)
        tk.Label(hdr, text=" Universal Cluster Optimizer",
                 fg=C["fg2"], bg=C["bg1"], font=FUI).pack(side="left")
        rh = tk.Frame(hdr, bg=C["bg1"])
        rh.pack(side="right", padx=12)
        mk_btn(rh, "Install MACE", self._install_mace,
               bg=C["bg3"], fg=C["fg2"], padx=10, pady=4).pack(side="right", padx=(6,0))
        self.lbl_mace = tk.Label(rh, text="MACE: checking…",
                                 fg=C["yel"], bg=C["bg1"], font=FUI)
        self.lbl_mace.pack(side="right")
        tk.Frame(self, bg=C["bg4"], height=1).pack(fill="x")

        # Preset toolbar
        ptb = tk.Frame(self, bg=C["bg"], pady=4)
        ptb.pack(fill="x", padx=8)
        tk.Label(ptb, text="Profile:", fg=C["fg3"], bg=C["bg"],
                 font=("Segoe UI",8)).pack(side="left", padx=(4,6))
        for name, col in [("Fast",C["grn"]),("Standard",C["blue"]),
                           ("Thorough",C["ora"]),("MACE+AL",C["pur"])]:
            mk_btn(ptb, name, lambda n=name: self._apply_preset(n),
                   bg=C["bg2"], fg=col,
                   font=("Segoe UI",8), padx=10, pady=2).pack(side="left", padx=2)
        tk.Label(ptb, text="  Ctrl+S save  Ctrl+O load  Ctrl+H history  F5 run  Esc stop",
                 fg=C["fg3"], bg=C["bg"],
                 font=("Segoe UI",8)).pack(side="right", padx=8)
        tk.Frame(self, bg=C["bg4"], height=1).pack(fill="x")
        body = tk.Frame(self, bg=C["bg"]); body.pack(fill="both", expand=True)
        self._build_sidebar(body)
        tk.Frame(body, bg=C["bg4"], width=1).pack(side="left", fill="y")
        self._build_content(body)

        # Status bar
        tk.Frame(self, bg=C["bg4"], height=1).pack(fill="x")
        bar = tk.Frame(self, bg=C["bg1"], height=52)
        bar.pack(fill="x", side="bottom"); bar.pack_propagate(False)
        self.btn_run = tk.Button(bar, text="▶  Run", command=self._run,
                                 bg=C["blue2"], fg="white",
                                 font=("Segoe UI",10,"bold"), relief="flat", bd=0,
                                 cursor="hand2", highlightthickness=0,
                                 activebackground="#1d4ed8", padx=24, pady=9)
        self.btn_run.pack(side="left", padx=(12,4), pady=9)
        self.btn_stop = tk.Button(bar, text="■  Stop", command=self._stop,
                                  bg=C["red"], fg="white", font=FBLD,
                                  relief="flat", bd=0, cursor="hand2",
                                  highlightthickness=0,
                                  activebackground="#c53030", padx=16,
                                  pady=9, state="disabled")
        self.btn_stop.pack(side="left", pady=9)
        self.lbl_status = tk.Label(bar, text="Ready", fg=C["grn"],
                                   bg=C["bg1"], font=FBLD)
        self.lbl_status.pack(side="left", padx=16)
        self.lbl_bin = tk.Label(bar, text="", fg=C["fg3"],
                                bg=C["bg1"], font=("Segoe UI",8))
        self.lbl_bin.pack(side="right", padx=12)

    # ── Settings persistence ──────────────────────────────────────────────────
    def _config_dict(self):
        """Serialise all tkinter variables to a plain dict."""
        try:
            geo = self.geometry()
        except Exception:
            geo = "1200x800"
        return {
            "formula":   self.vf.get(),   "pot":     self.vp.get(),
            "runs":      self.vruns.get(),"gens":    self.vgens.get(),
            "workers":   self.vwrk.get(), "seed":    self.vseed.get(),
            "bh_temp":   self.vbht.get(), "bh_decay":self.vbhd.get(),
            "step":      self.vstp.get(), "nfrac":   self.vnfr.get(),
            "graphprob": self.vgrp.get(), "bondcut": self.vbnd.get(),
            "tablam":    self.vtlm.get(), "tabdel":  self.vtld.get(),
            "gp":        self.bgp.get(),  "kappa":   self.vkap.get(),
            "cands":     self.vcnd.get(), "ard":     self.bard.get(),
            "ngd":       self.bngd.get(), "pop":     self.bpop.get(),
            "popsize":   self.vpop.get(), "homotop": self.bhtop.get(),
            "eigenvec":  self.beig.get(), "rmsd":    self.vrms.get(),
            "maxiso":    self.vmxi.get(), "ewin":    self.vewn.get(),
            "pt":        self.bpt.get(),  "pt_rep":  self.vptrp.get(),
            "pt_tmin":   self.vptlo.get(),"pt_tmax": self.vpthi.get(),
            "active":    self.bact.get(), "oracle":  self.vorc.get(),
            "outdir":    self.vout.get(), "verbose": self.bvrb.get(),
            "geometry":  geo,
            "history":   self._run_history[-20:],
        }

    def _apply_config(self, d):
        """Apply a config dict to all variables (ignores unknown keys)."""
        def s(var, key):
            if key in d: var.set(d[key])
        def b(var, key):
            if key in d: var.set(bool(d[key]))
        s(self.vf,"formula");   s(self.vp,"pot")
        s(self.vruns,"runs");   s(self.vgens,"gens")
        s(self.vwrk,"workers"); s(self.vseed,"seed")
        s(self.vbht,"bh_temp"); s(self.vbhd,"bh_decay")
        s(self.vstp,"step");    s(self.vnfr,"nfrac")
        s(self.vgrp,"graphprob");s(self.vbnd,"bondcut")
        s(self.vtlm,"tablam");  s(self.vtld,"tabdel")
        b(self.bgp,"gp");       s(self.vkap,"kappa")
        s(self.vcnd,"cands");   b(self.bard,"ard")
        b(self.bngd,"ngd");     b(self.bpop,"pop")
        s(self.vpop,"popsize"); b(self.bhtop,"homotop")
        b(self.beig,"eigenvec");s(self.vrms,"rmsd")
        s(self.vmxi,"maxiso");  s(self.vewn,"ewin")
        b(self.bpt,"pt");       s(self.vptrp,"pt_rep")
        s(self.vptlo,"pt_tmin");s(self.vpthi,"pt_tmax")
        b(self.bact,"active");  s(self.vorc,"oracle")
        s(self.vout,"outdir");  b(self.bvrb,"verbose")
        if "history" in d:
            self._run_history = d["history"]
        if "geometry" in d:
            try:
                self.geometry(d["geometry"])
            except Exception:
                pass

    def _save_config(self):
        try:
            CONFIG_FILE.write_text(json.dumps(self._config_dict(), indent=2))
        except Exception:
            pass   # non-fatal

    def _load_config(self):
        try:
            if CONFIG_FILE.exists():
                self._apply_config(json.loads(CONFIG_FILE.read_text()))
        except Exception:
            pass   # non-fatal: use defaults

    def _save_config_as(self):
        path = filedialog.asksaveasfilename(
            title="Save configuration",
            defaultextension=".json",
            initialfile="atlas_config.json",
            filetypes=[("JSON config","*.json"),("All","*.*")])
        if path:
            try:
                Path(path).write_text(json.dumps(self._config_dict(), indent=2))
                self._log_("ok", f"Config saved: {path}")
            except Exception as e:
                messagebox.showerror("Save error", str(e))

    def _load_config_from(self):
        path = filedialog.askopenfilename(
            title="Load configuration",
            filetypes=[("JSON config","*.json"),("All","*.*")])
        if path:
            try:
                self._apply_config(json.loads(Path(path).read_text()))
                self._log_("ok", f"Config loaded: {path}")
            except Exception as e:
                messagebox.showerror("Load error", str(e))

    # ── Preset profiles ────────────────────────────────────────────────────────
    PRESETS = {
        "Fast": {
            "runs":"4", "gens":"150", "workers":"4",
            "gp":True,  "ard":True,  "ngd":True,
            "pop":True, "eigenvec":True, "pt":False, "active":False,
            "kappa":"2.0", "cands":"3", "popsize":"10",
        },
        "Standard": {
            "runs":"8", "gens":"300", "workers":"4",
            "gp":True,  "ard":True,  "ngd":True,
            "pop":True, "eigenvec":True, "pt":False, "active":False,
            "kappa":"2.0", "cands":"5", "popsize":"20",
        },
        "Thorough": {
            "runs":"16","gens":"600", "workers":"4",
            "gp":True,  "ard":True,  "ngd":True,
            "pop":True, "eigenvec":True, "pt":True,
            "pt_rep":"4","pt_tmin":"0.3","pt_tmax":"5.0",
            "active":False, "kappa":"2.5", "cands":"8", "popsize":"30",
        },
        "MACE+AL": {
            "runs":"8", "gens":"200", "workers":"4",
            "pot":"MACE",
            "gp":True,  "ard":True,  "ngd":True,
            "pop":True, "eigenvec":True, "pt":False,
            "active":True, "oracle":"100",
            "kappa":"2.0", "cands":"5", "popsize":"20",
        },
    }

    def _apply_preset(self, name):
        self._apply_config(self.PRESETS[name])
        self._log_("info", f"Profile loaded: {name}")

    # ── Input validation ───────────────────────────────────────────────────────
    PARAM_RULES = [
        # (var, name, type, lo, hi)
        ("vruns",  "Runs",          int,   1,    256),
        ("vgens",  "Generations",   int,   10,   100000),
        ("vwrk",   "Workers",       int,   1,    256),
        ("vseed",  "Seed",          int,   0,    2**32),
        ("vbht",   "BH temperature",float, 0.01, 100.0),
        ("vbhd",   "BH decay",      float, 0.90, 1.0),
        ("vstp",   "Step fraction", float, 0.01, 5.0),
        ("vnfr",   "Atom fraction", float, 0.01, 1.0),
        ("vgrp",   "Graph-move prob",float,0.0,  1.0),
        ("vbnd",   "Bond cutoff",   float, 0.5,  5.0),
        ("vtlm",   "Taboo lambda",  float, 0.0,  1.0),
        ("vtld",   "Taboo delta",   float, 0.01, 5.0),
        ("vkap",   "GP kappa",      float, 0.1,  20.0),
        ("vcnd",   "GP candidates", int,   1,    100),
        ("vpop",   "Pop size",      int,   2,    200),
        ("vrms",   "RMSD thresh",   float, 0.01, 5.0),
        ("vmxi",   "Max isomers",   int,   1,    500),
        ("vewn",   "Energy window", float, 0.01, 10.0),
        ("vptrp",  "PT replicas",   int,   2,    16),
        ("vptlo",  "PT T min",      float, 0.01, 100.0),
        ("vpthi",  "PT T max",      float, 0.01, 1000.0),
        ("vorc",   "Oracle budget", int,   1,    10000),
    ]

    def _clear_field_errors(self):
        """Reset all entry borders to normal."""
        for e in self._entry_reg.values():
            try:
                e.config(highlightbackground=C["bg4"], highlightcolor=C["blue"],
                         highlightthickness=1)
            except Exception:
                pass

    def _flag_field(self, var):
        """Highlight a specific entry widget in red."""
        e = self._entry_reg.get(id(var))
        if e:
            try:
                e.config(highlightbackground=C["red"], highlightcolor=C["red"],
                         highlightthickness=2)
            except Exception:
                pass

    def _validate(self):
        """Validate all numeric fields.
        Highlights invalid fields in red inline, then shows a summary dialog.
        Returns True only if all fields pass.
        """
        self._clear_field_errors()
        errors = []
        for attr, name, typ, lo, hi in self.PARAM_RULES:
            var = getattr(self, attr)
            try:
                v = typ(var.get())
                if not (lo <= v <= hi):
                    errors.append(f"{name}: {v} outside [{lo}, {hi}]")
                    self._flag_field(var)
            except ValueError:
                errors.append(f"{name}: '{var.get()}' is not a valid {typ.__name__}")
                self._flag_field(var)

        if not self.vf.get().strip():
            errors.append("Formula: cannot be empty")

        try:
            if self.bpt.get() and float(self.vptlo.get()) >= float(self.vpthi.get()):
                errors.append("PT T min must be less than PT T max")
                self._flag_field(self.vptlo)
                self._flag_field(self.vpthi)
        except ValueError:
            pass

        if errors:
            messagebox.showerror(
                "Invalid parameters",
                "Please fix the highlighted fields:\n\n• " +
                "\n• ".join(errors))
            return False

        self._clear_field_errors()
        return True

    # ── ETA updater ────────────────────────────────────────────────────────────
    def _update_eta(self):
        """Called from _handle_json progress; updates status bar with ETA."""
        if not self._run_start_time or not self.best_per_run:
            return
        elapsed = time.time() - self._run_start_time
        # Average fraction complete across all active runs
        fracs = []
        for run, pts in self.best_per_run.items():
            tot = self.run_total.get(run, self.total_gens)
            if pts:
                fracs.append(pts[-1][0] / tot)
        if not fracs:
            return
        avg_frac = sum(fracs) / len(fracs)
        if avg_frac > 0.02:
            eta_total = elapsed / avg_frac
            eta_remain = max(0, eta_total - elapsed)
            m, s = divmod(int(eta_remain), 60)
            eta_str = f"{m}m {s:02d}s" if m else f"{s}s"
            self.lbl_status.config(
                text=f"● Running  {elapsed:.0f}s elapsed  ETA ≈ {eta_str}",
                fg=C["yel"])

    # ── Run history ────────────────────────────────────────────────────────────
    def _history_window(self):
        """Show a small window with the last N runs."""
        win = tk.Toplevel(self)
        win.title("Run history")
        win.configure(bg=C["bg"])
        win.geometry("620x320")
        win.resizable(True, True)

        tk.Label(win, text="Recent runs", fg=C["blue"], bg=C["bg"],
                 font=FH1).pack(anchor="w", padx=12, pady=(10,4))

        cols = ("Formula", "Potential", "E_min", "Isomers", "Time")
        tree_frame = tk.Frame(win, bg=C["bg"])
        tree_frame.pack(fill="both", expand=True, padx=12, pady=(0,8))

        sty = ttk.Style()
        sty.configure("H.Treeview", background=C["bg2"], fieldbackground=C["bg2"],
                       foreground=C["fg"], rowheight=22, font=FMONO)
        sty.configure("H.Treeview.Heading", background=C["bg3"],
                       foreground=C["fg2"], font=FUI)
        sty.map("H.Treeview", background=[("selected", C["blue2"])])

        tree = ttk.Treeview(tree_frame, columns=cols, show="headings",
                             style="H.Treeview")
        widths = [100, 80, 140, 80, 80]
        for col, w in zip(cols, widths):
            tree.heading(col, text=col)
            tree.column(col, width=w, anchor="center")

        sc = tk.Scrollbar(tree_frame, command=tree.yview)
        tree.configure(yscrollcommand=sc.set)
        tree.pack(side="left", fill="both", expand=True)
        sc.pack(side="right", fill="y")

        for rec in reversed(self._run_history):
            tree.insert("", "end", values=(
                rec.get("formula","?"),
                rec.get("pot","?"),
                f"{rec.get('best_e', 0.0):.8f}",
                rec.get("n_iso","?"),
                f"{rec.get('elapsed', 0.0):.1f}s",
            ))

        def _reload(evt):
            sel = tree.selection()
            if not sel:
                return
            vals = tree.item(sel[0])["values"]
            self.vf.set(str(vals[0]))
            self.vp.set(str(vals[1]))
            win.destroy()

        tree.bind("<Double-1>", _reload)
        tk.Label(win, text="Double-click a row to reload formula + potential",
                 fg=C["fg3"], bg=C["bg"], font=("Segoe UI",8)).pack(pady=(0,6))

    # ── Keyboard shortcuts ─────────────────────────────────────────────────────
    def _bind_shortcuts(self):
        self.bind("<F5>",        lambda _: self._run())
        self.bind("<Escape>",    lambda _: self._stop() if self.running else None)
        self.bind("<Control-s>", lambda _: self._save_config_as())
        self.bind("<Control-o>", lambda _: self._load_config_from())
        self.bind("<Control-h>", lambda _: self._history_window())
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _on_close(self):
        self._save_config()
        self.destroy()

    # ── Build sidebar ─────────────────────────────────────────────────────────
    def _build_sidebar(self, parent):
        sb = tk.Frame(parent, bg=C["bg1"], width=310)
        sb.pack(side="left", fill="y"); sb.pack_propagate(False)

        cv = tk.Canvas(sb, bg=C["bg1"], highlightthickness=0, bd=0, width=306)
        cv.pack(side="left", fill="both", expand=True)
        sc = tk.Scrollbar(sb, orient="vertical", command=cv.yview)
        sc.pack(side="right", fill="y")
        cv.configure(yscrollcommand=sc.set)
        inner = tk.Frame(cv, bg=C["bg1"])
        wid = cv.create_window((0,0), window=inner, anchor="nw")
        inner.bind("<Configure>", lambda e: cv.configure(scrollregion=cv.bbox("all")))
        cv.bind("<Configure>", lambda e: cv.itemconfig(wid, width=e.width))
        cv.bind_all("<MouseWheel>",
                    lambda e: cv.yview_scroll(int(-1*(e.delta/120)),"units"))

        P = {"padx":10, "pady":3}

        # ── Cluster ─────────────────────────────────────────────────
        divider(inner,"Cluster").pack(fill="x",**P)
        ff = tk.Frame(inner, bg=C["bg1"]); ff.pack(fill="x",padx=10,pady=4)
        mk_lbl(ff,"Formula").pack(side="left")
        self.fe = tk.Entry(ff, textvariable=self.vf, bg=C["bg3"], fg=C["blue"],
                           insertbackground=C["fg"], relief="flat", bd=0,
                           font=("Consolas",13,"bold"), width=13,
                           highlightthickness=1, highlightbackground=C["bg4"],
                           highlightcolor=C["blue"])
        self.fe.pack(side="left", padx=6)
        tip(self.fe, HELP["formula"])
        self.fe.bind("<Return>", lambda _: self._run())

        pf = tk.Frame(inner, bg=C["bg1"]); pf.pack(fill="x",padx=10,pady=2)
        mk_lbl(pf,"Potential").pack(side="left")
        sty = ttk.Style()
        sty.configure("P.TCombobox", fieldbackground=C["bg3"],
                      background=C["bg3"], foreground=C["fg"],
                      selectbackground=C["bg3"], selectforeground=C["blue"])
        cb = ttk.Combobox(pf, textvariable=self.vp, style="P.TCombobox",
                          values=["AUTO","LJ","MORSE","GUPTA","UFF","MACE"],
                          font=FUI, width=8, state="readonly")
        cb.pack(side="left", padx=6)
        tip(cb, HELP["potential"])
        mk_lbl(pf,"AUTO = best for formula", fg=C["fg3"],
               font=("Segoe UI",8)).pack(side="left")

        mk_lbl(inner,"Quick presets:", fg=C["fg3"],
               font=("Segoe UI",8)).pack(anchor="w",padx=10,pady=(4,2))
        qf = tk.Frame(inner, bg=C["bg1"]); qf.pack(fill="x",padx=10,pady=(0,6))
        for i,(fm,pt) in enumerate([("Ar13","LJ"),("LJ38","LJ"),
                                     ("Ni13","GUPTA"),("Au55","GUPTA"),
                                     ("Au32Ag6","GUPTA"),("Ni12Cu1","GUPTA"),
                                     ("Fe12O18","UFF"),("Pt13","MACE")]):
            mk_btn(qf, fm, lambda f=fm,p=pt:(self.vf.set(f),self.vp.set(p)),
                   bg=C["bg2"], fg=C["blue"],
                   font=("Segoe UI",8), padx=6, pady=2
                   ).grid(row=i//4, column=i%4, padx=2, pady=2, sticky="ew")

        # ── Search ──────────────────────────────────────────────────
        divider(inner,"Search").pack(fill="x",**P)
        for lbl,var,tipkey,hint in [
            ("Runs",        self.vruns,"runs",  "Independent restarts"),
            ("Generations", self.vgens,"gens",  "Steps per run"),
            ("Workers",     self.vwrk, "workers","Parallel CPU threads"),
            ("Random seed", self.vseed,"seed",  "0 = random"),
        ]:
            r = param_row(inner,lbl,var,tip=hint,registry=self._entry_reg)
            r.pack(fill="x",padx=10,pady=2)
            tip(r, HELP[tipkey])

        # ── Basin-hopping ────────────────────────────────────────────
        divider(inner,"Basin-hopping").pack(fill="x",**P)
        for lbl,var,tipkey,hint in [
            ("BH temperature", self.vbht,"bht",    "Metropolis accept temp"),
            ("Cooling rate",   self.vbhd,"bhd",    "--bh-decay per step"),
            ("Step fraction",  self.vstp,"step",   "Perturbation size / σ"),
            ("Atoms per step", self.vnfr,"nfrac",  "Fraction displaced"),
            ("Graph-move prob",self.vgrp,"grprob", "Bond-topology moves"),
            ("Bond cutoff",    self.vbnd,"bondcut","×σ, defines connectivity"),
        ]:
            r = param_row(inner,lbl,var,tip=hint,registry=self._entry_reg)
            r.pack(fill="x",padx=10,pady=2)
            tip(r, HELP[tipkey])

        # ── Taboo memory ─────────────────────────────────────────────
        divider(inner,"Taboo memory").pack(fill="x",**P)
        for lbl,var,tipkey,hint in [
            ("Lambda",self.vtlm,"tablam","Repulsion strength"),
            ("Delta", self.vtld,"tabdel","Repulsion width"),
        ]:
            r = param_row(inner,lbl,var,tip=hint,registry=self._entry_reg)
            r.pack(fill="x",padx=10,pady=2)
            tip(r, HELP[tipkey])

        # ── GP surrogate ─────────────────────────────────────────────
        divider(inner,"GP surrogate").pack(fill="x",**P)
        chk_gp = mk_chk(inner,"Use GP surrogate",self.bgp)
        chk_gp.pack(anchor="w",padx=10,pady=2)
        tip(chk_gp, "Enable the Gaussian Process Bayesian surrogate.\n"
                    "Scores random candidate moves and directs search\n"
                    "toward low-energy unexplored regions (LCB acquisition).")
        for lbl,var,tipkey,hint in [
            ("GP κ (kappa)",  self.vkap,"kappa",  "Explore/exploit balance"),
            ("GP candidates", self.vcnd,"cands",  "Moves scored per step"),
            ("Pop size",      self.vpop,"popsize","Elite crossover pool"),
        ]:
            r = param_row(inner,lbl,var,tip=hint,registry=self._entry_reg)
            r.pack(fill="x",padx=10,pady=2)
            tip(r, HELP[tipkey])

        # ── Isomer filter ────────────────────────────────────────────
        divider(inner,"Isomer filter").pack(fill="x",**P)
        for lbl,var,tipkey,hint in [
            ("RMSD threshold",self.vrms,"rmsd",  "Similarity cutoff ×σ"),
            ("Max isomers",   self.vmxi,"maxiso","Isomers to keep"),
            ("Energy window", self.vewn,"ewin",  "Fraction above min"),
        ]:
            r = param_row(inner,lbl,var,tip=hint,registry=self._entry_reg)
            r.pack(fill="x",padx=10,pady=2)
            tip(r, HELP[tipkey])

        # ── Algorithms ───────────────────────────────────────────────
        divider(inner,"Algorithms").pack(fill="x",**P)
        for var,lbl,tipkey in [
            (self.bard, "ARD adaptive descriptor (S4)",   "ard"),
            (self.bngd, "Natural gradient hyperparams (S4)","ngd"),
            (self.bpop, "Population crossover (S4)",      "pop"),
            (self.bhtop,"Homotop search (heteroatomic)",  "homotop"),
            (self.beig, "Lanczos eigenvector escape (S3)","eigenvec"),
        ]:
            c = mk_chk(inner,lbl,var)
            c.pack(anchor="w",padx=10,pady=2)
            tip(c, HELP[tipkey])

        # ── Parallel tempering ───────────────────────────────────────
        divider(inner,"Parallel tempering").pack(fill="x",**P)
        c_pt = mk_chk(inner,"Enable",self.bpt)
        c_pt.pack(anchor="w",padx=10,pady=2)
        tip(c_pt, HELP["pt"])
        ptg = tk.Frame(inner,bg=C["bg1"]); ptg.pack(fill="x",padx=10,pady=2)
        for col,(lbl,var,tipkey) in enumerate([
            ("Replicas",self.vptrp,"ptreplicas"),
            ("T min",   self.vptlo,"ptmin"),
            ("T max",   self.vpthi,"ptmax"),
        ]):
            lw = mk_lbl(ptg,lbl,font=("Segoe UI",8))
            lw.grid(row=0,column=col*2,padx=(0,2))
            ew = mk_entry(ptg,var,w=5)
            ew.grid(row=0,column=col*2+1,padx=(0,8))
            tip(ew, HELP[tipkey])

        # ── Active learning ──────────────────────────────────────────
        divider(inner,"Active learning (MACE oracle)").pack(fill="x",**P)
        c_al = mk_chk(inner,"Enable (requires MACE)",self.bact)
        c_al.pack(anchor="w",padx=10,pady=2)
        tip(c_al, HELP["active"])
        r_orc = param_row(inner,"Oracle budget",self.vorc,tip="max MACE calls",registry=self._entry_reg)
        r_orc.pack(fill="x",padx=10,pady=2)
        tip(r_orc, HELP["oracle"])

        # ── Output ───────────────────────────────────────────────────
        divider(inner,"Output").pack(fill="x",**P)
        od = tk.Frame(inner,bg=C["bg1"]); od.pack(fill="x",padx=10,pady=2)
        mk_lbl(od,"Output dir",width=14,anchor="w").pack(side="left")
        oe = mk_entry(od,self.vout,w=15); oe.pack(side="left",padx=4)
        tip(oe, HELP["outdir"])
        mk_btn(od,"…",self._browse,padx=4,pady=1).pack(side="left")
        mk_chk(inner,"Verbose statistics (--verbose)",self.bvrb).pack(
            anchor="w",padx=10,pady=(2,14))



    # ── Content area ──────────────────────────────────────────────────────────
    def _build_content(self, parent):
        area = tk.Frame(parent, bg=C["bg"]); area.pack(fill="both",expand=True)
        sty = ttk.Style()
        sty.configure("M.TNotebook", background=C["bg"], borderwidth=0)
        sty.configure("M.TNotebook.Tab", background=C["bg1"], foreground=C["fg2"],
                      font=FUI, padding=[16,7])
        sty.map("M.TNotebook.Tab",
                background=[("selected",C["bg"])],
                foreground=[("selected",C["blue"])])
        nb = ttk.Notebook(area, style="M.TNotebook")
        nb.pack(fill="both",expand=True)
        self.nb = nb
        self.t_mon = tk.Frame(nb,bg=C["bg"])
        self.t_res = tk.Frame(nb,bg=C["bg"])
        self.t_log = tk.Frame(nb,bg=C["bg"])
        nb.add(self.t_mon, text="  Monitor  ")
        nb.add(self.t_res, text="  Results  ")
        nb.add(self.t_log, text="  Log  ")
        self._build_monitor()
        self._build_results()
        self._build_log()

    # ── Monitor tab ───────────────────────────────────────────────────────────
    def _build_monitor(self):
        p = self.t_mon
        # Stats row
        sr = tk.Frame(p,bg=C["bg"]); sr.pack(fill="x",padx=12,pady=(10,6))
        self._sc = {}
        for key,lbl,col in [("formula","Cluster",C["blue"]),
                              ("best","Global min",C["grn"]),
                              ("delta","Improvement",C["yel"]),
                              ("basins","New basins",C["pur"]),
                              ("gen","Generation",C["ora"]),
                              ("nruns","Runs",C["fg2"])]:
            cd = tk.Frame(sr, bg=C["bg2"], highlightthickness=1,
                          highlightbackground=C["bg4"])
            cd.pack(side="left", padx=3)
            tk.Label(cd,text=lbl,fg=C["fg2"],bg=C["bg2"],
                     font=("Segoe UI",8)).pack(padx=12,pady=(5,0))
            lb = tk.Label(cd,text="—",fg=col,bg=C["bg2"],
                          font=("Consolas",11,"bold"),width=14)
            lb.pack(padx=12,pady=(0,5))
            self._sc[key] = lb
        # Progress bar
        pf = tk.Frame(p,bg=C["bg"]); pf.pack(fill="x",padx=12,pady=(0,4))
        sty = ttk.Style()
        sty.configure("A.Horizontal.TProgressbar", troughcolor=C["bg2"],
                      background=C["blue"], bordercolor=C["bg2"])
        self.pbar = ttk.Progressbar(pf, style="A.Horizontal.TProgressbar",
                                    orient="horizontal", mode="determinate",
                                    maximum=100)
        self.pbar.pack(fill="x")
        # Canvas
        self.cv = tk.Canvas(p, bg=C["bg2"], highlightthickness=0)
        self.cv.pack(fill="both",expand=True,padx=12,pady=(0,12))
        self.cv.bind("<Configure>", lambda _: self.after(10,self._redraw))

    def _sc_set(self, key, val):
        if key in self._sc: self._sc[key].config(text=str(val))

    # ── Results tab ───────────────────────────────────────────────────────────
    def _build_results(self):
        p = self.t_res
        p.columnconfigure(1, weight=1); p.rowconfigure(0, weight=1)

        # List
        L = tk.Frame(p,bg=C["bg"]); L.grid(row=0,column=0,sticky="nsew",padx=(12,0),pady=10)
        L.rowconfigure(1,weight=1)
        tk.Label(L,text="Isomers",fg=C["blue"],bg=C["bg"],font=FH1).grid(row=0,column=0,sticky="w",pady=(0,6))
        lo = tk.Frame(L,bg=C["bg2"],highlightthickness=1,highlightbackground=C["bg4"])
        lo.grid(row=1,column=0,sticky="nsew"); lo.rowconfigure(0,weight=1); lo.columnconfigure(0,weight=1)
        sc = tk.Scrollbar(lo); sc.grid(row=0,column=1,sticky="ns")
        self.iso_lb = tk.Listbox(lo,bg=C["bg2"],fg=C["fg"],font=FMONO,
                                  selectbackground=C["blue2"],selectforeground="white",
                                  yscrollcommand=sc.set,relief="flat",
                                  activestyle="none",width=30,bd=0,highlightthickness=0)
        self.iso_lb.grid(row=0,column=0,sticky="nsew")
        sc.config(command=self.iso_lb.yview)
        self.iso_lb.bind("<<ListboxSelect>>",self._show_iso)
        br = tk.Frame(L,bg=C["bg"]); br.grid(row=2,column=0,sticky="ew",pady=(6,0))
        mk_btn(br,"Export XYZ",self._export,padx=10,pady=4).pack(side="left")
        mk_btn(br,"Open folder",self._open_dir,padx=10,pady=4).pack(side="left",padx=4)

        # Separator
        tk.Frame(p,bg=C["bg4"],width=1).grid(row=0,column=0,sticky="ns",padx=(262,0))

        # Detail
        R = tk.Frame(p,bg=C["bg"]); R.grid(row=0,column=1,sticky="nsew",padx=10,pady=10)
        R.rowconfigure(1,weight=1); R.columnconfigure(0,weight=1)
        self.iso_hdr = tk.Label(R,text="Select an isomer",fg=C["fg2"],
                                bg=C["bg"],font=FBLD,anchor="w")
        self.iso_hdr.grid(row=0,column=0,sticky="ew",pady=(0,6))
        to = tk.Frame(R,bg=C["bg2"],highlightthickness=1,highlightbackground=C["bg4"])
        to.grid(row=1,column=0,sticky="nsew")
        to.rowconfigure(0,weight=1); to.columnconfigure(0,weight=1)
        sc2 = tk.Scrollbar(to); sc2.grid(row=0,column=1,sticky="ns")
        self.iso_tx = tk.Text(to,bg=C["bg2"],fg=C["fg"],font=FMONO,
                               yscrollcommand=sc2.set,relief="flat",
                               state="disabled",wrap="none",bd=0,
                               highlightthickness=0,padx=10,pady=8)
        self.iso_tx.grid(row=0,column=0,sticky="nsew")
        sc2.config(command=self.iso_tx.yview)
        self.iso_tx.tag_config("hdr",foreground=C["blue"],font=FBLD)
        self.iso_tx.tag_config("num",foreground=C["fg2"],font=FMONO)
        self.iso_tx.tag_config("dim",foreground=C["fg3"],font=FMONO)

    # ── Log tab ───────────────────────────────────────────────────────────────
    def _build_log(self):
        p = self.t_log
        p.rowconfigure(0,weight=1); p.columnconfigure(0,weight=1)
        sc = tk.Scrollbar(p); sc.grid(row=0,column=1,sticky="ns",padx=(0,4),pady=8)
        self.log = tk.Text(p,bg=C["bg2"],fg=C["fg2"],font=FMONO,
                           yscrollcommand=sc.set,relief="flat",state="disabled",
                           wrap="word",bd=0,highlightthickness=0,padx=12,pady=8)
        self.log.grid(row=0,column=0,sticky="nsew",padx=(12,0),pady=8)
        sc.config(command=self.log.yview)
        for tag,col in [("ok",C["grn"]),("info",C["blue"]),("warn",C["yel"]),
                         ("error",C["red"]),("result",C["pur"]),("cmd",C["fg3"]),
                         ("dim",C["fg3"])]:
            self.log.tag_config(tag,foreground=col)
        btn_row = tk.Frame(p, bg=C["bg"]); btn_row.grid(row=1,column=0,sticky="e",padx=12,pady=(0,8))
        mk_btn(btn_row,"Copy log",self._copy_log,bg=C["bg3"],fg=C["fg2"],padx=10,pady=3).pack(side="left",padx=(0,4))
        mk_btn(btn_row,"Clear",   self._clrlog,  bg=C["bg3"],fg=C["fg2"],padx=10,pady=3).pack(side="left")

    # ── Command builder (every flag) ──────────────────────────────────────────
    def _build_cmd(self):
        cmd = [self.binary, self.vf.get(),
               "--pot",           self.vp.get(),
               "--runs",          self.vruns.get(),
               "--gens",          self.vgens.get(),
               "--workers",       self.vwrk.get(),
               "--seed",          self.vseed.get(),
               "--bh-temp",       self.vbht.get(),
               "--bh-decay",      self.vbhd.get(),
               "--stress-step",   self.vstp.get(),
               "--stress-nfrac",  self.vnfr.get(),
               "--graph-prob",    self.vgrp.get(),
               "--bond-cutoff",   self.vbnd.get(),
               "--taboo-lambda",  self.vtlm.get(),
               "--taboo-delta",   self.vtld.get(),
               "--gp-kappa",      self.vkap.get(),
               "--gp-cands",      self.vcnd.get(),
               "--pop-size",      self.vpop.get(),
               "--rmsd-thresh",   self.vrms.get(),
               "--max-isomers",   self.vmxi.get(),
               "--energy-window", self.vewn.get(),
               "--output-dir",    self.vout.get(),
               "--mace-script",   self.mace_script,
               "--progress"]
        if not self.bgp.get():   cmd.append("--no-gp")
        if not self.bard.get():  cmd.append("--no-ard")
        if not self.bngd.get():  cmd.append("--no-ngd")
        if not self.bpop.get():  cmd.append("--no-pop")
        if not self.bhtop.get(): cmd.append("--no-homotop")
        if not self.beig.get():  cmd.append("--no-eigenvec")
        if self.bvrb.get():      cmd.append("--verbose")
        if self.bpt.get():
            cmd += ["--pt","--pt-replicas",self.vptrp.get(),
                    "--pt-tmin",self.vptlo.get(),"--pt-tmax",self.vpthi.get()]
        if self.bact.get():
            cmd += ["--active","--oracle-budget",self.vorc.get()]
        return cmd

    # ── Run / Stop ────────────────────────────────────────────────────────────
    def _run(self):
        if self.running: return
        if not os.path.exists(self.binary):
            messagebox.showerror("Binary not found",
                f"atlas.exe not found:\n{self.binary}\n\nRun ATLAS.bat first.")
            return
        if not self._validate():     # ← validation gate
            return
        self._run_start_time = time.time()
        self.running = True
        self.isomers = []; self.best_per_run = {}
        self.run_total = {}; self.global_min = None
        self.iso_lb.delete(0,"end"); self.pbar["value"] = 0
        self._sc_set("formula",self.vf.get())
        for k in ("best","delta","basins","gen","nruns"):
            self._sc_set(k,"—")
        self.cv.delete("all")
        self.btn_run.config(state="disabled"); self.btn_stop.config(state="normal")
        self.lbl_status.config(text="● Running",fg=C["yel"])
        self.nb.select(self.t_mon)
        cmd = self._build_cmd()
        self._log_("cmd","$ "+" ".join(cmd))
        def worker():
            try:
                self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT,
                                             text=True, bufsize=1)
                for raw in self.proc.stdout:
                    line = raw.rstrip()
                    if not line: continue
                    if line.startswith("{"):
                        self.after(0,self._handle_json,line)
                    else:
                        self.after(0,self._logline,line)
                self.proc.wait()
            except Exception as ex:
                self.after(0,self._log_,"error",str(ex))
            finally:
                self.after(0,self._done)
        threading.Thread(target=worker,daemon=True).start()

    def _stop(self):
        if self.proc:
            try: self.proc.terminate()
            except: pass
        self._done()

    def _done(self):
        self.running = False
        self._run_start_time = None
        self.btn_run.config(state="normal"); self.btn_stop.config(state="disabled")
        self.lbl_status.config(text="Ready",fg=C["grn"])
        self._save_config()     # persist settings on every run completion

    # ── JSON handler ──────────────────────────────────────────────────────────
    def _handle_json(self, line):
        try: d = json.loads(line)
        except: return
        t = d.get("type","")

        if t == "header":
            self.n_runs = d.get("runs",8)
            self.total_gens = d.get("gens",300)
            fm = d.get("formula","?"); N = d.get("N","?"); pot = d.get("potential","?")
            self._sc_set("formula",f"{fm}  N={N}")
            self._log_("info",f"Start: {fm} N={N} pot={pot} "
                              f"{self.n_runs}×{self.total_gens} gens")

        elif t == "progress":
            run = d.get("run",0); gen = d.get("gen",0)
            tot = d.get("total",self.total_gens)
            best = d.get("best",0.0); bas = d.get("basins",0)
            if run not in self.best_per_run: self.best_per_run[run] = []
            self.best_per_run[run].append((gen,best))
            self.run_total[run] = tot
            gb = min(pts[-1][1] for pts in self.best_per_run.values() if pts)
            if self.global_min is None or gb < self.global_min:
                self.global_min = gb
            pts = self.best_per_run[run]
            dstr = f"{pts[-1][1]-pts[0][1]:+.5f}" if len(pts)>=2 else "—"
            fracs = [pts[-1][0]/self.run_total.get(r,tot)
                     for r,pts in self.best_per_run.items() if pts]
            self.pbar["value"] = (sum(fracs)/len(fracs)*100) if fracs else 0
            self._sc_set("best",f"{gb:.8f}")
            self._sc_set("delta",dstr)
            self._sc_set("basins",str(bas))
            self._sc_set("gen",f"{gen}/{tot}")
            self._sc_set("nruns",f"{len(self.best_per_run)}/{self.n_runs}")
            self._update_eta()          # ← ETA in status bar
            self._redraw()

        elif t == "result":
            be = d.get("best_energy",0.0); ni = d.get("n_isomers",0)
            el = d.get("elapsed",0.0); pot = d.get("potential","?")
            self._sc_set("best",f"{be:.8f}"); self.pbar["value"] = 100
            self._log_("result",f"DONE  E_min={be:.8f}  {ni} isomers  "
                                f"{el:.1f}s  pot={pot}")
            # Record in history
            self._run_history.append({
                "formula": self.vf.get(), "pot": pot,
                "best_e": be, "n_iso": ni, "elapsed": el,
            })
            self._load_iso(); self.nb.select(self.t_res)

    def _logline(self, line):
        u = line.upper()
        if "[OK]" in u or "[ATLAS]" in u: self._log_("ok",line)
        elif "ERROR" in u: self._log_("error",line)
        elif "WARN" in u: self._log_("warn",line)
        elif line.startswith("==="): self._log_("result",line)
        else: self._log_("info",line)

    def _log_(self, lvl, msg):
        tag = {"ok":"ok","info":"info","warn":"warn","error":"error",
               "result":"result","cmd":"cmd"}.get(lvl,"info")
        ts  = time.strftime("%H:%M:%S")
        self.log.config(state="normal")
        self.log.insert("end", f"{ts} ", "dim")
        self.log.insert("end", f"[{lvl.upper():6s}] {msg}\n", tag)
        self.log.see("end"); self.log.config(state="disabled")

    def _clrlog(self):
        self.log.config(state="normal"); self.log.delete("1.0","end")
        self.log.config(state="disabled")

    def _copy_log(self):
        text = self.log.get("1.0", "end").strip()
        self.clipboard_clear()
        self.clipboard_append(text)
        self._log_("info", "Log copied to clipboard")

    # ── Convergence plot ──────────────────────────────────────────────────────
    def _redraw(self):
        cv = self.cv; cv.delete("all")
        W = cv.winfo_width(); H = cv.winfo_height()
        if W<80 or H<80: return
        cv.create_rectangle(0,0,W,H,fill=C["bg2"],outline="")
        PL,PR,PT,PB = 80,20,28,40
        if not self.best_per_run:
            cv.create_text(W//2,H//2,text="Waiting for data…",
                           fill=C["fg3"],font=FUI); return
        all_e = [e for pts in self.best_per_run.values() for _,e in pts]
        max_g = max(self.run_total.get(r,pts[-1][0])
                    for r,pts in self.best_per_run.items() if pts) or 1
        elo,ehi = min(all_e),max(all_e)
        if ehi==elo: ehi=elo+abs(elo)*0.01 or 1.0
        def px(g): return PL+(g/max_g)*(W-PL-PR)
        def py(e): return PT+(1-(e-elo)/(ehi-elo))*(H-PT-PB)
        # Grid
        for i in range(5):
            frac=i/4; y=py(elo+frac*(ehi-elo))
            cv.create_line(PL,y,W-PR,y,fill=C["bg4"],dash=(3,8))
            cv.create_text(PL-5,y,text=f"{elo+frac*(ehi-elo):.4f}",
                           fill=C["fg3"],font=("Consolas",7),anchor="e")
            x=px(frac*max_g)
            cv.create_line(x,PT,x,H-PB,fill=C["bg4"],dash=(3,8))
            cv.create_text(x,H-PB+6,text=str(int(frac*max_g)),
                           fill=C["fg3"],font=("Consolas",7))
        cv.create_line(PL,PT,PL,H-PB,fill=C["bg4"],width=1)
        cv.create_line(PL,H-PB,W-PR,H-PB,fill=C["bg4"],width=1)
        # Per-run staircase traces
        for rid,pts in sorted(self.best_per_run.items()):
            if len(pts)<2: continue
            col=TRACES[rid%len(TRACES)]
            coords=[]
            for k,(g,e) in enumerate(pts):
                if k==0: coords.extend([px(g),py(e)])
                else:
                    prev_y=py(pts[k-1][1])
                    coords.extend([px(g),prev_y,px(g),py(e)])
            if len(coords)>=4: cv.create_line(coords,fill=col,width=1.5)
        # Global min envelope
        all_pts = sorted([(g,e) for pts in self.best_per_run.values()
                          for g,e in pts],key=lambda x:x[0])
        if len(all_pts)>=2:
            rmin=all_pts[0][1]; gc=[px(all_pts[0][0]),py(rmin)]
            for g,e in all_pts[1:]:
                if e<rmin:
                    gc.extend([px(g),py(rmin),px(g),py(e)]); rmin=e
                else:
                    gc.extend([px(g),py(rmin)])
            if len(gc)>=4: cv.create_line(gc,fill=C["grn"],width=2.5)
        # Title + legend
        cv.create_text(PL+6,12,text="Energy convergence  —  all runs",
                       fill=C["fg2"],font=FUI,anchor="w")
        lx=PL+6
        for rid in sorted(self.best_per_run.keys())[:8]:
            cv.create_line(lx,H-10,lx+16,H-10,fill=TRACES[rid%len(TRACES)],width=1.5)
            cv.create_text(lx+20,H-10,text=f"r{rid+1}",fill=C["fg3"],
                           font=("Consolas",7),anchor="w"); lx+=44
        cv.create_line(lx,H-10,lx+16,H-10,fill=C["grn"],width=2.5)
        cv.create_text(lx+20,H-10,text="global min",fill=C["grn"],
                       font=("Consolas",7,"bold"),anchor="w")

    # ── Isomer list ───────────────────────────────────────────────────────────
    def _load_iso(self):
        path = os.path.join(self.vout.get(),f"{self.vf.get()}_isomers.xyz")
        if not os.path.exists(path):
            self._log_("warn",f"Not found: {path}"); return
        self.isomers=[]
        try:
            with open(path) as f:
                while True:
                    nl=f.readline()
                    if not nl: break
                    nl=nl.strip()
                    if not nl: continue
                    N=int(nl); com=f.readline().strip()
                    coords=[f.readline().split() for _ in range(N)]
                    self.isomers.append({"N":N,"comment":com,"coords":coords})
        except Exception as ex:
            self._log_("error",f"Parse: {ex}"); return
        energies=[self._eparse(iso["comment"]) for iso in self.isomers]
        e0=energies[0] if energies and isinstance(energies[0],float) else None
        self.iso_lb.delete(0,"end")
        for i,(iso,e) in enumerate(zip(self.isomers,energies)):
            rank="★" if i==0 else f"#{i+1:3d}"
            es=f"{e:.8f}" if isinstance(e,float) else "  ?"
            ds="" if (i==0 or not isinstance(e,float) or e0 is None) else f"  Δ{e-e0:+.5f}"
            self.iso_lb.insert("end",f" {rank}  {es}{ds}")
        if self.isomers:
            self.iso_lb.selection_set(0); self._show_iso(None)
        self._log_("ok",f"Loaded {len(self.isomers)} isomers")

    def _eparse(self, comment):
        for tok in comment.split():
            if tok.startswith("Energy="):
                try: return float(tok[7:])
                except: pass
        return "?"

    def _show_iso(self, _):
        sel=self.iso_lb.curselection()
        if not sel: return
        idx=sel[0]
        if idx>=len(self.isomers): return
        iso=self.isomers[idx]; e=self._eparse(iso["comment"])
        es=f"{e:.10f}" if isinstance(e,float) else "?"
        label="Global minimum ★" if idx==0 else f"Isomer #{idx+1}"
        self.iso_hdr.config(text=f"{label}  —  E = {es}",
                             fg=C["grn"] if idx==0 else C["blue"])
        t=self.iso_tx; t.config(state="normal"); t.delete("1.0","end")
        t.insert("end",f"{'#':>5s}  {'El':^4s}  {'X (Å)':>14s}  {'Y (Å)':>14s}  {'Z (Å)':>14s}\n","hdr")
        t.insert("end","─"*62+"\n","hdr")
        for k,parts in enumerate(iso["coords"]):
            if len(parts)<4: continue
            elem=parts[0]
            try: x,y,z=float(parts[1]),float(parts[2]),float(parts[3])
            except: continue
            col=ecol(elem)
            tag_e=f"el_{elem}"
            if tag_e not in t.tag_names(): t.tag_config(tag_e,foreground=col,font=FMONO)
            t.insert("end",f"{k+1:5d}  ","num")
            t.insert("end",f"{elem:^4s}",tag_e)
            t.insert("end",f"  {x:14.8f}  {y:14.8f}  {z:14.8f}\n","num")
        t.insert("end",f"\nComment: {iso['comment']}\n","dim")
        t.insert("end",f"Atoms:   {iso['N']}\n","dim")
        t.config(state="disabled"); t.see("1.0")

    def _export(self):
        if not self.isomers:
            messagebox.showinfo("No isomers","Run a search first."); return
        src=os.path.join(self.vout.get(),f"{self.vf.get()}_isomers.xyz")
        if not os.path.exists(src):
            messagebox.showerror("Not found",f"File not found:\n{src}"); return
        dst=filedialog.asksaveasfilename(title="Export isomers",
                                          defaultextension=".xyz",
                                          initialfile=f"{self.vf.get()}_isomers.xyz",
                                          filetypes=[("XYZ","*.xyz"),("All","*.*")])
        if dst: shutil.copy(src,dst); self._log_("ok",f"Exported: {dst}")

    def _open_dir(self):
        d=self.vout.get(); os.makedirs(d,exist_ok=True)
        if sys.platform=="win32": os.startfile(d)
        elif sys.platform=="darwin": subprocess.Popen(["open",d])
        else: subprocess.Popen(["xdg-open",d])

    def _browse(self):
        d=filedialog.askdirectory(title="Output directory")
        if d: self.vout.set(d)

    # ── MACE installer ────────────────────────────────────────────────────────
    def _install_mace(self):
        script=os.path.join(os.path.dirname(self.mace_script),"install_mace.py")
        def _do():
            self.lbl_mace.config(text="MACE: installing…",fg=C["yel"])
            self._log_("info","Installing MACE-MP-0 (5–15 min)…")
            if not os.path.exists(script):
                self._log_("error",f"install_mace.py not found: {script}")
                self.lbl_mace.config(text="MACE: error",fg=C["red"])
                return
            try:
                pr=subprocess.Popen([sys.executable,script],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT,
                                    text=True,bufsize=1)
                for raw in pr.stdout:
                    raw=raw.strip()
                    if not raw: continue
                    try:
                        d2=json.loads(raw)
                        st=d2.get("status",""); msg=d2.get("message",raw)
                        # If install_mace.py found a compatible Python, remember it
                        if "python" in d2 and d2["python"]:
                            self.mace_python = d2["python"]
                            self._log_("info",f"MACE will use: {self.mace_python}")
                        if st in("success","already_installed"):
                            self.lbl_mace.config(text="MACE: installed ✓",fg=C["grn"])
                            self._log_("ok",msg)
                        elif st=="error":
                            self.lbl_mace.config(text="MACE: error",fg=C["red"])
                            self._log_("error",msg)
                        else:
                            self._log_("info",msg)
                    except:
                        self._log_("info",raw)
                pr.wait()
            except Exception as ex:
                self.lbl_mace.config(text="MACE: failed",fg=C["red"])
                self._log_("error",str(ex))
            threading.Thread(target=self._chk_mace,daemon=True).start()
        threading.Thread(target=_do,daemon=True).start()

    def _chk_mace(self):
        """Check MACE in current Python, then try py -3.11 / py -3.12 fallback."""
        # First try current Python
        try:
            r=subprocess.run([sys.executable,"-c","import mace;print(mace.__version__)"],
                             capture_output=True,text=True,timeout=15)
            if r.returncode==0:
                self.lbl_mace.config(text=f"MACE {r.stdout.strip()} ✓",fg=C["grn"])
                return
        except Exception:
            pass
        # Fallback: check py -3.11 and py -3.12
        import shutil
        if shutil.which("py"):
            for ver in ("3.11","3.12","3.10"):
                try:
                    r=subprocess.run(["py",f"-{ver}","-c",
                                      "import mace;print(mace.__version__)"],
                                     capture_output=True,text=True,timeout=15)
                    if r.returncode==0:
                        v=r.stdout.strip()
                        self.lbl_mace.config(text=f"MACE {v} (py {ver}) ✓",fg=C["grn"])
                        # Remember which Python has MACE
                        r2=subprocess.run(["py",f"-{ver}","-c",
                                           "import sys;print(sys.executable)"],
                                          capture_output=True,text=True,timeout=10)
                        if r2.returncode==0:
                            self.mace_python=r2.stdout.strip()
                        return
                except Exception:
                    continue
        self.lbl_mace.config(text="MACE: not installed",fg=C["yel"])

    def _show_shortcuts(self):
        win = tk.Toplevel(self)
        win.title("Keyboard shortcuts")
        win.configure(bg=C["bg"])
        win.geometry("340x220")
        win.resizable(False,False)
        shortcuts = [
            ("F5",        "Run search"),
            ("Escape",    "Stop running search"),
            ("Ctrl+S",    "Save configuration to file"),
            ("Ctrl+O",    "Load configuration from file"),
            ("Ctrl+H",    "Show run history"),
            ("Enter",     "Run (when formula field focused)"),
            ("Hover",     "Show parameter help tooltip"),
        ]
        tk.Label(win, text="Keyboard shortcuts", fg=C["blue"], bg=C["bg"],
                 font=FH1).pack(anchor="w", padx=14, pady=(12,8))
        for key, desc in shortcuts:
            r = tk.Frame(win, bg=C["bg"])
            r.pack(fill="x", padx=14, pady=1)
            tk.Label(r, text=key, fg=C["yel"], bg=C["bg"], font=FMONO,
                     width=10, anchor="w").pack(side="left")
            tk.Label(r, text=desc, fg=C["fg2"], bg=C["bg"],
                     font=FUI, anchor="w").pack(side="left")

    def _startup(self):
        if os.path.exists(self.binary):
            self.lbl_bin.config(text=f"atlas.exe: {self.binary}",fg=C["fg3"])
            self._log_("ok",f"Binary: {self.binary}")
        else:
            self.lbl_bin.config(text="atlas.exe: NOT FOUND",fg=C["red"])
            self._log_("warn",f"atlas.exe not found at {self.binary} — run ATLAS.bat first")
        threading.Thread(target=self._chk_mace,daemon=True).start()
        # Refresh history label if any runs already done
        if self._run_history:
            self._log_("info", f"Loaded {len(self._run_history)} previous runs from history")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="ATLAS-GP GUI")
    ap.add_argument("--binary", default="build/atlas.exe")
    ap.add_argument("--mace",   default="ml/mace_server.py")
    args = ap.parse_args()
    root = Path(__file__).parent.parent
    binary = args.binary if os.path.isabs(args.binary) else str(root/args.binary)
    mace   = args.mace   if os.path.isabs(args.mace)   else str(root/args.mace)
    ATLASApp(binary, mace).mainloop()

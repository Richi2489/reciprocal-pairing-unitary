"""
run_mu_scaling.py — Estudio de escalamiento de μ(d) para n≥4 (Ruta 1).
Ver docs/mu_scaling_prereg.md (reglas pre-registradas).

Malla n∈[4,5,6,8,10] × familias {vm, const(0.4), const_fixedS(Σ=2.0, control C2),
linear, random×3}. num_starts = min(max(200,10n²),400). Emite CSV incremental y
un reporte con gráficas de escala y análisis del confounder Σd_i.
"""

import csv
import os
import sys
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, ValueError):  # pragma: no cover
    pass

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from qgraphs.mu_general import (
    min_mu_unitary, von_mangoldt_diag, constant_diag, constant_diag_fixed_sum,
    linear_diag, random_diag,
)

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(RAIZ, "data", "mu_scaling_results.csv")
FIGDIR = os.path.join(RAIZ, "docs", "figures", "qgraphs")
REPORTE = os.path.join(RAIZ, "docs", "mu_scaling_report.md")

EPS0, DELTA = 0.05, 1e-3
N_LIST = [4, 5, 6, 8, 10]


def num_starts(n):
    return min(max(200, 10 * n * n), 400)


def familias(n):
    """(etiqueta, seed_reportado, d) para cada n."""
    out = [
        ("vm", 0, von_mangoldt_diag(n)),
        ("const", 0, constant_diag(n, 0.4)),
        ("const_fixedS", 0, constant_diag_fixed_sum(n, 2.0)),   # control C2
        ("linear", 0, linear_diag(n)),
    ]
    for s in (0, 1, 2):
        out.append(("random", s, random_diag(n, s)))
    return out


CAMPOS = ["n", "family", "seed", "sum_d", "mu", "p5", "median",
          "err_constraint", "n_starts", "n_feasible", "status"]


def correr():
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    with open(CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CAMPOS)
        w.writeheader()
        for n in N_LIST:
            ns = num_starts(n)
            for fam, seed, d in familias(n):
                r = min_mu_unitary(d, num_starts=ns, seed=0)
                if r["infeasible_apriori"]:
                    status, mu, p5, med, err, nf = "infeasible_apriori", "", "", "", "", 0
                elif not r["feasible"]:
                    status, mu, p5, med, err, nf = "infeasible_numeric", "", "", "", "", 0
                else:
                    mins = np.array(r["mins"])
                    mu = r["mu"]; p5 = float(np.percentile(mins, 5))
                    med = float(np.median(mins)); err = r["err"]
                    nf = r["n_starts_feasible"]
                    status = "ok"
                    if mins.min() < DELTA:
                        status = "FALSIFIED"          # algún mínimo bajo δ
                    # persistir el minimizador para autopsias post-hoc (lección del
                    # FALSIFIED de const_fixedS n=8, irrecuperable sin re-correr)
                    udir = os.path.join(RAIZ, "data", "mu_scaling_U")
                    os.makedirs(udir, exist_ok=True)
                    np.savez_compressed(
                        os.path.join(udir, f"n{n}_{fam}_s{seed}.npz"),
                        U=r["U"], d=d, mu=mu, err=err)
                row = dict(n=n, family=fam, seed=seed, sum_d=round(float(d.sum()), 4),
                           mu=mu, p5=p5, median=med, err_constraint=err,
                           n_starts=ns, n_feasible=nf, status=status)
                w.writerow(row)
                f.flush()
                print(f"n={n:2d} {fam:13s} seed={seed} Σd={d.sum():6.3f} "
                      f"mu={mu if mu=='' else round(mu,4)} status={status}", flush=True)


# ---------------------------------------------------------------------------
# Análisis + reporte
# ---------------------------------------------------------------------------


def _ajuste(ns, ys):
    """Ajuste log-lineal a potencia (y=a n^-alpha) y exponencial (y=a e^-bn).
    Devuelve (alpha, R2_pot, beta, R2_exp)."""
    ns = np.asarray(ns, float); ys = np.asarray(ys, float)
    if len(ns) < 2 or np.any(ys <= 0):
        return None
    ly = np.log(ys)
    # potencia: ln y = ln a - alpha ln n
    A = np.polyfit(np.log(ns), ly, 1); alpha = -A[0]
    r2p = _r2(ly, np.polyval(A, np.log(ns)))
    # exponencial: ln y = ln a - beta n
    B = np.polyfit(ns, ly, 1); beta = -B[0]
    r2e = _r2(ly, np.polyval(B, ns))
    return alpha, r2p, beta, r2e


def _r2(y, yhat):
    ss = np.sum((y - np.mean(y)) ** 2)
    return float(1 - np.sum((y - yhat) ** 2) / ss) if ss > 0 else 1.0


def reporte():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = list(csv.DictReader(open(CSV)))
    os.makedirs(FIGDIR, exist_ok=True)

    def rho_de(n, sum_d, mu):
        """Fraccion de emparejamiento forzado ρ = 2μ/(n−Σd) ∈ [0,1]. Ver mu_ruta2_budget.md."""
        P = n - sum_d
        return float(2.0 * mu / P) if P > 1e-9 else float("nan")

    def fiable(r):
        """Celda fiable si ≥50% de los arranques fueron factibles (n=8,10 no lo son)."""
        return int(r["n_feasible"]) >= 0.5 * int(r["n_starts"])

    def serie(fam, solo_fiables=False):
        pts = [(int(r["n"]), float(r["median"]), float(r["sum_d"]), float(r["mu"]),
                fiable(r))
               for r in rows if r["family"] == fam and r["status"] in ("ok", "FALSIFIED")
               and (fiable(r) or not solo_fiables)]
        pts.sort()
        return pts

    def serie_rho(fam, solo_fiables=False):
        """(n, ρ, fiable) usando μ del mejor mínimo, mediana sobre seeds si hay varias."""
        d = {}
        for r in rows:
            if r["family"] == fam and r["status"] in ("ok", "FALSIFIED") \
                    and (fiable(r) or not solo_fiables):
                n = int(r["n"])
                d.setdefault(n, []).append(
                    (rho_de(n, float(r["sum_d"]), float(r["mu"])), fiable(r)))
        return sorted((n, float(np.median([v for v, _ in lst])),
                       all(f for _, f in lst)) for n, lst in d.items())

    def _marca_no_fiable(ax_o_plt):
        """Sombrea la región n>6 (factibilidad colapsada, ver reporte §1/§4)."""
        ax_o_plt.axvspan(6.5, 10.5, color="gray", alpha=0.12,
                         label="n>6: NO fiable (factib.<50%)")

    fams_lisas = ["vm", "const", "const_fixedS", "linear"]

    # --- Fig 1: μ (mediana) vs n; línea sólo por puntos fiables, no-fiables huecos ---
    plt.figure(figsize=(7, 5))
    for fam in fams_lisas:
        p = serie(fam)
        if p:
            fi = [(a, b) for a, b, _, _, ok in p if ok]
            no = [(a, b) for a, b, _, _, ok in p if not ok]
            ln, = plt.plot([a for a, _ in fi], [b for _, b in fi], "o-", label=fam)
            if no:
                plt.plot([a for a, _ in no], [b for _, b in no], "o", mfc="none",
                         color=ln.get_color())
    _marca_no_fiable(plt)
    plt.axhline(EPS0, ls="--", c="gray", lw=0.8, label=f"ε₀={EPS0}")
    plt.xlabel("n"); plt.ylabel("μ (mediana)"); plt.title("μ(d) vs dimensión n")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "mu_scaling_vs_n.png"), dpi=130); plt.close()

    # --- Fig 2: confounder μ vs Σd (fiables llenos, no-fiables huecos) ---
    plt.figure(figsize=(7, 5))
    for fam in sorted(set(r["family"] for r in rows)):
        cel = [r for r in rows if r["family"] == fam
               and r["status"] in ("ok", "FALSIFIED")]
        fi = [(float(r["sum_d"]), float(r["mu"])) for r in cel if fiable(r)]
        no = [(float(r["sum_d"]), float(r["mu"])) for r in cel if not fiable(r)]
        if fi:
            sc = plt.scatter([a for a, _ in fi], [b for _, b in fi], label=fam, s=30)
            if no:
                plt.scatter([a for a, _ in no], [b for _, b in no], s=30,
                            facecolors="none", edgecolors=sc.get_facecolor())
        elif no:
            plt.scatter([a for a, _ in no], [b for _, b in no], s=30,
                        facecolors="none", label=f"{fam} (no fiable)")
    plt.xlabel("Σd_i"); plt.ylabel("μ")
    plt.title("Confounder: μ vs masa diagonal total (huecos = no fiable)")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "mu_vs_sumd.png"), dpi=130); plt.close()

    # --- Fig 3: control C2 — μ vs n con Σd fijo=2.0 (solo puntos fiables en línea) ---
    plt.figure(figsize=(7, 5))
    p = serie("const_fixedS")
    if p:
        fi = [(a, b) for a, b, _, _, ok in p if ok]
        no = [(a, b) for a, b, _, _, ok in p if not ok]
        plt.plot([a for a, _ in fi], [b for _, b in fi], "s-", c="crimson",
                 label="const, Σd=2.0 (C2)")
        if no:
            plt.plot([a for a, _ in no], [b for _, b in no], "s", mfc="none",
                     c="crimson")
    _marca_no_fiable(plt)
    plt.axhline(EPS0, ls="--", c="gray", lw=0.8)
    plt.xlabel("n"); plt.ylabel("μ (mediana)")
    plt.title("Control C2: μ vs n con Σd_i fijo (ojo: las entradas 2/n encogen)")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "mu_C2_fixedS.png"), dpi=130); plt.close()

    # --- Fig 4: ρ = 2μ/(n−Σd), fracción de emparejamiento forzado, vs n ---
    plt.figure(figsize=(7, 5))
    for fam in fams_lisas:
        pr = serie_rho(fam)
        if pr:
            fi = [(a, b) for a, b, ok in pr if ok]
            no = [(a, b) for a, b, ok in pr if not ok]
            ln, = plt.plot([a for a, _ in fi], [b for _, b in fi], "o-", label=fam)
            if no:
                plt.plot([a for a, _ in no], [b for _, b in no], "o", mfc="none",
                         color=ln.get_color())
    _marca_no_fiable(plt)
    plt.ylim(0, 1); plt.xlabel("n")
    plt.ylabel("ρ = 2μ/(n−Σd)  (fracción de emparejamiento forzado)")
    plt.title("Obstrucción normalizada ρ vs n  (¿piso positivo o →0?)")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "rho_vs_n.png"), dpi=130); plt.close()

    # --- ajustes de escala ---
    ajustes = {}
    for fam in fams_lisas:
        p = serie(fam)
        if len(p) >= 2:
            ajustes[fam] = _ajuste([a for a, *_ in p], [b for _, b, *_ in p])

    # --- redacción ---
    L = []
    L.append("# Estudio de escalamiento de μ(d) para n ≥ 4 — reporte\n")
    L.append("_open_math_lab (banco Hilbert–Pólya). Ruta 1 (numérica). "
             "Reglas: `docs/mu_scaling_prereg.md`. Generado por `scripts/run_mu_scaling.py`._\n")
    L.append("_Evidencia numérica (cota superior de μ + positividad multi-start), "
             "NO prueba. La prueba es la Ruta 2 (O1/O2 analíticas)._\n")

    L.append("\n## Resultados (μ = mejor mínimo factible; mediana entre arranques)\n")
    L.append("ρ = 2μ/(n−Σd) = fracción del presupuesto fuera-de-diagonal forzada a emparejarse "
             "(ver `docs/mu_ruta2_budget.md`; observable adimensional de la obstrucción).\n")
    L.append("| n | familia | seed | Σd | μ | ρ | p5 | mediana | err | #feas | estado |")
    L.append("|---|---------|------|-----|-----|-----|-----|---------|-----|-------|--------|")
    for r in rows:
        mu = "" if r["mu"] == "" else f"{float(r['mu']):.4f}"
        rho = "" if r["mu"] == "" else f"{rho_de(int(r['n']), float(r['sum_d']), float(r['mu'])):.3f}"
        p5 = "" if r["p5"] == "" else f"{float(r['p5']):.4f}"
        md = "" if r["median"] == "" else f"{float(r['median']):.4f}"
        er = "" if r["err_constraint"] == "" else f"{float(r['err_constraint']):.1e}"
        L.append(f"| {r['n']} | {r['family']} | {r['seed']} | {r['sum_d']} | {mu} | {rho} "
                 f"| {p5} | {md} | {er} | {r['n_feasible']} | {r['status']} |")

    L.append("\n## Ajustes de escala (mediana de μ vs n)\n")
    L.append("| familia | potencia μ~n^(−α): α | R² | exponencial μ~e^(−βn): β | R² |")
    L.append("|---------|----------------------|-----|--------------------------|-----|")
    for fam, a in ajustes.items():
        if a:
            alpha, r2p, beta, r2e = a
            L.append(f"| {fam} | {alpha:.3f} | {r2p:.3f} | {beta:.3f} | {r2e:.3f} |")

    # --- veredictos pre-registrados ---
    ok = [r for r in rows if r["status"] in ("ok", "FALSIFIED")]
    falsif = [r for r in rows if r["status"] == "FALSIFIED"]
    p5_min = min((float(r["p5"]) for r in ok), default=float("nan"))
    infeas = [r for r in rows if "infeasible" in r["status"]]
    etiquetas_infeas = sorted(set(
        "{}(n={},seed={})".format(r["family"], r["n"], r["seed"]) for r in infeas))
    lista_infeas = ", ".join(etiquetas_infeas) if etiquetas_infeas else "—"
    L.append("\n## Veredictos (reglas pre-registradas)\n")
    L.append(f"- **Positividad** (p5 > ε₀={EPS0}): p5 mínimo sobre todas las celdas factibles "
             f"= {p5_min:.4f} → {'SÍ, evidencia μ(d)>0' if p5_min > EPS0 else 'NO en alguna celda'}.")
    L.append(f"- **Falsación** (algún mínimo < δ={DELTA}): "
             f"{'DISPARADA en ' + str(len(falsif)) + ' celda(s)' if falsif else 'no disparada (0 celdas)'}.")
    L.append(f"- **Infactibilidad detectada**: {len(infeas)} celda(s) ({lista_infeas}).")

    # --- confounder / C2 ---
    # --- ρ(n): el observable normalizado (pregunta central de Ruta 2) ---
    L.append("\n## Obstrucción normalizada ρ(n) = 2μ/(n−Σd)\n")
    L.append("Cota superior DEMOSTRADA: μ(d) ≤ (n−Σd)/2 ⇒ ρ∈[0,1] (`docs/mu_ruta2_budget.md`). "
             "μ(d)>0 ⟺ ρ>0. Pregunta: ¿ρ tiene piso positivo universal o →0 con n?\n")
    L.append("| familia | ρ(n) sobre n creciente | ¿tendencia? |")
    L.append("|---------|------------------------|-------------|")
    for fam in fams_lisas:
        pr = serie_rho(fam, solo_fiables=True)   # la tendencia SOLO sobre celdas fiables
        if len(pr) >= 2:
            serie_txt = ", ".join(f"n={a}:{b:.3f}" for a, b, *_ in pr)
            tend = "DECRECE" if pr[-1][1] < pr[0][1] - 0.02 else (
                "satura/plano" if abs(pr[-1][1] - pr[0][1]) <= 0.02 else "crece")
            L.append(f"| {fam} | {serie_txt} | {tend} |")
    L.append("\nFig. `rho_vs_n.png`. Si ρ→0 la obstrucción-por-lazo se DILUYE en dimensión alta "
             "(consecuencia directa para si un grafo finito puede realizar las amplitudes de ζ); "
             "si satura en un piso positivo, es candidata a cota universal (objetivo O1/O2).\n")

    L.append("\n## Confounder Σd_i y control C2\n")
    c2 = serie("const_fixedS")
    if len(c2) >= 2:
        mus = [b for _, b, *_ in c2]
        rango = max(mus) - min(mus)
        L.append(f"- Con **Σd_i fijo = 2.0** (const_fixedS), μ va de {min(mus):.4f} a {max(mus):.4f} "
                 f"sobre n∈{[a for a,*_ in c2]} (rango {rango:.4f}). "
                 f"{'μ CAMBIA con n a Σd fijo → hay efecto DIMENSIONAL puro.' if rango > 0.02 else 'μ ~estable a Σd fijo → el escalamiento aparente es sobre todo de MAGNITUD (Σd), no de dimensión.'}")
    L.append("- Fig. `mu_vs_sumd.png`: dispersión μ vs Σd (la correlación cruda que vimos en n=4).")
    L.append("- Fig. `mu_C2_fixedS.png`: μ vs n con Σd fijo (dimensión aislada) — la brújula para la Ruta 2.")

    L.append("\n## Figuras\n- `docs/figures/qgraphs/mu_scaling_vs_n.png`\n"
             "- `docs/figures/qgraphs/mu_vs_sumd.png`\n"
             "- `docs/figures/qgraphs/mu_C2_fixedS.png`\n")

    open(REPORTE, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print(f"\nReporte: {REPORTE}", flush=True)


if __name__ == "__main__":
    _desconocidos = [a for a in sys.argv[1:] if a not in ("--figures-only", "--report-only")]
    if _desconocidos:
        sys.exit(f"argumento(s) desconocido(s): {_desconocidos}. Uso: python -m scripts.run_mu_scaling "
                 "[--report-only | --figures-only]. SIN argumentos re-corre TODA la malla y reescribe "
                 "data/mu_scaling_results.csv y docs/mu_scaling_report.md.")
    if "--figures-only" in sys.argv:
        # regenera SOLO las figuras; NO reescribe docs/mu_scaling_report.md, que es
        # el reporte final editado a mano (5 secciones) y no debe aplastarse
        _rep = REPORTE
        REPORTE = os.path.join(RAIZ, "data", "_reporte_auto_descartable.md")
        reporte()
        os.remove(REPORTE)
        REPORTE = _rep
    elif "--report-only" in sys.argv:
        reporte()   # regenera reporte+figuras desde el CSV existente, sin re-correr la malla
    else:
        correr()
        reporte()

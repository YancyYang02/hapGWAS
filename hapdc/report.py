
"""hapdc.report -- tables and dependency-free SVG figures."""
import csv
import os

import numpy as np

__all__ = ["write_windows_table", "write_summary", "write_manhattan_svg", "write_qq_svg"]

COL = dict(A_joint="#1f6feb", P1_sidak="#2ea043", D_acat="#d1242f",
           HDCb="#8250df", HDCt="#bf8700")


def write_windows_table(path, res, trait):
    """One row per window: every statistic's observed p, plus window annotations."""
    P = res["P"]
    keys = [k for k in ("A_joint", "P1_sidak", "D_acat", "HDCb", "HDCt") if k in P]
    n = P[keys[0]].shape[0]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["trait", "window", "marker_start", "nk", "ncell"] + keys)
        for i in range(n):
            w.writerow([trait, i, res["starts"][i], res["nk"][i], res["ncell"][i]]
                       + ["%.6g" % P[k][i, 0] for k in keys])
    return path


def write_summary(path, rows):
    """rows: list of dict with a shared key set, written as CSV."""
    if not rows:
        return None
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in keys})
    return path


def _axes(w, h, ml, mr, mt, mb):
    return ml, h - mb, w - mr, mt


def write_qq_svg(path, pvals, label="", width=380, height=380):
    p = np.sort(np.asarray(pvals, float)[np.isfinite(pvals)])
    n = len(p)
    if n == 0:
        return None
    exp = -(np.log10((np.arange(1, n + 1) - 0.5) / n))
    obs = -np.log10(np.clip(p, 1e-300, 1))
    lo, hi = 0, max(exp.max(), obs.max()) * 1.05
    ml, bot, rt, top = _axes(width, height, 52, 14, 26, 44)
    X = lambda v: ml + (v - lo) / (hi - lo) * (rt - ml)
    Y = lambda v: bot - (v - lo) / (hi - lo) * (bot - top)
    s = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="Helvetica,Arial" font-size="11">' % (width, height)]
    s.append('<rect width="100%" height="100%" fill="#fff"/>')
    s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#ccc" stroke-dasharray="3 3"/>' % (X(lo), Y(lo), X(hi), Y(hi)))
    s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#333"/>' % (ml, bot, rt, bot))
    s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#333"/>' % (ml, bot, ml, top))
    pts = " ".join("%.1f,%.1f" % (X(a), Y(b)) for a, b in zip(exp, obs))
    s.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.6"/>' % (pts, COL.get(label, "#1f6feb")))
    s.append('<text x="%d" y="16" text-anchor="middle" font-size="12">%s</text>' % (width // 2, label))
    s.append('<text x="%d" y="%d" text-anchor="middle">expected -log10 p</text>' % (width // 2, height - 10))
    s.append('<text x="14" y="%d" text-anchor="middle" transform="rotate(-90 14 %d)">observed -log10 p</text>' % ((top + bot) // 2, (top + bot) // 2))
    s.append('</svg>')
    open(path, "w").write("\n".join(s))
    return path


def write_manhattan_svg(path, pvals, starts, chrom_of, label="", width=760, height=300):
    p = np.asarray(pvals, float)
    ok = np.isfinite(p) & (p > 0)
    if not ok.any():
        return None
    y = -np.log10(np.clip(p, 1e-300, 1))
    ml, bot, rt, top = _axes(width, height, 48, 14, 26, 40)
    Y = lambda v: bot - v / (y[ok].max() * 1.08) * (bot - top)
    base = np.median(np.sort(y[ok])[: max(10, ok.sum() // 2)])
    thr = -np.log10(0.05 / max(1, ok.sum()))
    s = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="Helvetica,Arial" font-size="11">' % (width, height)]
    s.append('<rect width="100%" height="100%" fill="#fff"/>')
    idx = np.where(ok)[0]
    for j, i in enumerate(idx):
        x = ml + (i + 0.5) / (len(p)) * (rt - ml)
        s.append('<rect x="%.1f" y="%.1f" width="1.6" height="%.1f" fill="%s"/>'
                 % (x, Y(y[i]), max(bot - Y(y[i]), 0.8), COL.get(label, "#1f6feb")))
    s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#999" stroke-dasharray="4 3"/>' % (ml, Y(thr), rt, Y(thr)))
    s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#333"/>' % (ml, bot, rt, bot))
    s.append('<text x="%d" y="16" text-anchor="middle" font-size="12">%s</text>' % (width // 2, label))
    s.append('<text x="14" y="%d" text-anchor="middle" transform="rotate(-90 14 %d)">-log10 p</text>' % ((top + bot) // 2, (top + bot) // 2))
    s.append('<text x="%d" y="%d" text-anchor="end" font-size="10" fill="#666">dashed: Bonferroni</text>' % (rt, top + 2))
    s.append('</svg>')
    open(path, "w").write("\n".join(s))
    return path

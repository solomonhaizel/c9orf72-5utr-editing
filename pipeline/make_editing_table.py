#!/usr/bin/env python3
"""
Build the canonical per-replicate editing table (editing_by_rep.csv).

This is the single source of truth for every plotted editing value. Each row is
one construct x site x dose x biological replicate, with the guide raw editing,
the matched same-run no-guide (UTP) baseline, the net (guide - no-guide), and
coverage at the target. R plotting scripts read this table; no editing value is
transcribed by hand.

Editing = 100 * G / (A + G) at the target adenosine.
Net      = guide_raw - matched_noguide_raw (same sequencing run, per replicate).

Inputs (relative to repo root; raw reads are in the SRA, see README):
  .var files from each run's Genotyper output, and the per-position CSVs for the
  runs processed before the .var pipeline existed (C9-Main, ACTB reps 1-2).

Target positions:
  A74  = ref 107 in the 253-nt amplicon (= biological 74)
  A88  = ref 121 (= biological 88)
  ACTB = position 113 in its amplicon (edited A of the premature TAG)
"""
import csv, os, sys

VAR = "data/raw_var"           # folder holding each run's .var files (see data/README)
CSV = "data/raw_perposition"   # folder holding the legacy per-position CSVs

def edit_cov_var(path, pos):
    for r in csv.DictReader(open(path), delimiter="\t"):
        if int(r["pos"]) == pos:
            A, G = int(r["A"]), int(r["G"])
            cov = int(r["reads_all"])
            return (100*G/(A+G) if A+G else None), cov
    return None, None

def edit_cov_csv(path, sample, pos, acol="A_count", gcol="G_count", covcol=None):
    d = {(x["sample"], int(x["position"])): x for x in csv.DictReader(open(path))}
    r = d.get((sample, pos))
    if not r: return None, None
    A, G = float(r[acol]), float(r[gcol])
    cov = int(float(r[covcol])) if covcol and covcol in r else int(A+G)
    return (100*G/(A+G) if A+G else None), cov

# manifest: one entry per construct/site/dose/rep. 'src' picks the reader.
# guide + noguide identify the same-run matched pair.
M = [
 # ---- A74 (target ref 107 in .var; biological 74 in C9-Main csv) ----
 dict(construct="A74", site="A74", dose="0.3 uM", rep=1, run="run1",   src="csv", path=f"{CSV}/C9-Main.csv",  guide="A743P",  noguide="UTP", pos=74),
 dict(construct="A74", site="A74", dose="1 uM",   rep=1, run="run1",   src="csv", path=f"{CSV}/C9-Main.csv",  guide="A741uP", noguide="UTP", pos=74),
 dict(construct="A74", site="A74", dose="0.3 uM", rep=2, run="run2",   src="var", dir=f"{VAR}/run2",  guide="Plate1-A2", noguide="Plate1-A4", pos=107),
 dict(construct="A74", site="A74", dose="1 uM",   rep=2, run="run2",   src="var", dir=f"{VAR}/run2",  guide="Plate1-A1", noguide="Plate1-A4", pos=107),
 dict(construct="A74", site="A74", dose="0.3 uM", rep=3, run="run3",   src="var", dir=f"{VAR}/run3",  guide="Plate1-A2", noguide="Plate1-A4", pos=107),
 dict(construct="A74", site="A74", dose="1 uM",   rep=3, run="run3",   src="var", dir=f"{VAR}/run3",  guide="Plate1-A1", noguide="Plate1-A4", pos=107),
 dict(construct="A74", site="A74", dose="0.3 uM", rep=4, run="run4",   src="var", dir=f"{VAR}/run4",  guide="Plate1-A2", noguide="Plate1-A3", pos=107),
 dict(construct="A74", site="A74", dose="1 uM",   rep=4, run="run4",   src="var", dir=f"{VAR}/run4",  guide="Plate1-A1", noguide="Plate1-A3", pos=107),
 # ---- A88 (target ref 121; biological 88) ----
 dict(construct="A88", site="A88", dose="0.3 uM", rep=1, run="run1",   src="csv", path=f"{CSV}/C9-Main.csv",  guide="A883P",  noguide="UTP", pos=88),
 dict(construct="A88", site="A88", dose="1 uM",   rep=1, run="run1",   src="csv", path=f"{CSV}/C9-Main.csv",  guide="A881uP", noguide="UTP", pos=88),
 dict(construct="A88", site="A88", dose="0.3 uM", rep=2, run="run6",  src="var", dir=f"{VAR}/run6", guide="Plate1-A2", noguide="Plate1-A7", pos=121),
 dict(construct="A88", site="A88", dose="1 uM",   rep=2, run="run6",  src="var", dir=f"{VAR}/run6", guide="Plate1-A1", noguide="Plate1-A7", pos=121),
 dict(construct="A88", site="A88", dose="0.3 uM", rep=3, run="run5",  src="var", dir=f"{VAR}/run5", guide="Plate1-A3", noguide="Plate1-A7", pos=121),
 dict(construct="A88", site="A88", dose="1 uM",   rep=3, run="run5",  src="var", dir=f"{VAR}/run5", guide="Plate1-A1", noguide="Plate1-A7", pos=121),
 dict(construct="A88", site="A88", dose="0.3 uM", rep=4, run="run5",  src="var", dir=f"{VAR}/run5", guide="Plate1-A4", noguide="Plate1-A8", pos=121),
 dict(construct="A88", site="A88", dose="1 uM",   rep=4, run="run5",  src="var", dir=f"{VAR}/run5", guide="Plate1-A2", noguide="Plate1-A8", pos=121),
 # ---- ACTB reps 1-2 (target pos 113; dose swap: A1uP=0.3uM, A3P=1uM) ----
 dict(construct="ACTB", site="ACTB", dose="0.3 uM", rep=1, run="run1", src="acsv", path=f"{CSV}/ACTB-rep1.csv", guide="A1uP", noguide="A0P", pos=113),
 dict(construct="ACTB", site="ACTB", dose="1 uM",   rep=1, run="run1", src="acsv", path=f"{CSV}/ACTB-rep1.csv", guide="A3P",  noguide="A0P", pos=113),
 dict(construct="ACTB", site="ACTB", dose="0.3 uM", rep=2, run="run1", src="acsv", path=f"{CSV}/ACTB-rep2.csv", guide="A1uP", noguide="A0P", pos=113),
 dict(construct="ACTB", site="ACTB", dose="1 uM",   rep=2, run="run1", src="acsv", path=f"{CSV}/ACTB-rep2.csv", guide="A3P",  noguide="A0P", pos=113),
 # ACTB rep3 from the run7 run (target position 113).
 dict(construct="ACTB", site="ACTB", dose="0.3 uM",   rep=3, run="run7", src="var", dir=f"{VAR}/run7", guide="Plate1-A6", noguide="Plate1-A2", pos=113),
 dict(construct="ACTB", site="ACTB", dose="1 uM",     rep=3, run="run7", src="var", dir=f"{VAR}/run7", guide="Plate1-A4", noguide="Plate1-A2", pos=113),
 dict(construct="ACTB", site="ACTB", dose="scramble", rep=3, run="run7", src="var", dir=f"{VAR}/run7", guide="Plate1-A8", noguide="Plate1-A2", pos=113),
]

def main():
    out = []
    for m in M:
        if m["src"] == "var":
            g, gc = edit_cov_var(f"{m['dir']}/{m['guide']}.var", m["pos"])
            b, bc = edit_cov_var(f"{m['dir']}/{m['noguide']}.var", m["pos"])
        elif m["src"] == "csv":
            g, gc = edit_cov_csv(m["path"], m["guide"], m["pos"], covcol="total_depth")
            b, bc = edit_cov_csv(m["path"], m["noguide"], m["pos"], covcol="total_depth")
        elif m["src"] == "acsv":
            g, gc = edit_cov_csv(m["path"], m["guide"], m["pos"], covcol="reads_all")
            b, bc = edit_cov_csv(m["path"], m["noguide"], m["pos"], covcol="reads_all")
        net = round(g - b, 4) if (g is not None and b is not None) else ""
        out.append([m["construct"], m["site"], m["dose"], m["rep"], m["run"],
                    round(g,4) if g is not None else "", gc,
                    round(b,4) if b is not None else "", bc, net])

    with open("data/editing_by_rep.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["construct","site","dose","rep","run",
                    "guide_raw_pct","guide_cov","noguide_raw_pct",
                    "noguide_cov","net_pct"])
        w.writerows(out)
    print(f"wrote data/editing_by_rep.csv ({len(out)} rows)")

if __name__ == "__main__":
    main()

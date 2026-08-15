#!/usr/bin/env python3
"""
A-to-I editing at the ACTB positive-control target from Genotyper .var output.

The ACTB reporter carries a premature TAG stop; ADAR editing of that adenosine
(A->G, restoring TGG read-through) is scored at a FIXED reference position,
default 113 -- ACTB has no ATGTGA motif, so the C9 script (analyze_c9_editing.py)
must NOT be used here. This script takes the target position explicitly and
validates it is a reference adenosine before scoring.

Editing = 100 * G / (A + G) at the target. No baseline is subtracted here; the
matched same-run no-guide (A0P) subtraction is done per replicate in
make_editing_table.py.

Usage:
  python analyze_actb_editing.py RUN_DIR [--run NAME] [--pos 113] [--min-cov 10000]
  RUN_DIR contains the Plate1-*.var files (ref0.fa optional).
"""
import argparse, csv, glob, os, sys

def counts_at(var_path, pos):
    for r in csv.DictReader(open(var_path), delimiter="\t"):
        if int(r["pos"]) == pos:
            return {b: int(r[b]) for b in "ACGT"} | {"reads_all": int(r["reads_all"]), "ref": r["ref"]}
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--run", default=None)
    ap.add_argument("--pos", type=int, default=113, help="ACTB target reference position")
    ap.add_argument("--min-cov", type=int, default=10000)
    args = ap.parse_args()
    run = args.run or os.path.basename(os.path.normpath(args.run_dir))

    vfiles = sorted(glob.glob(os.path.join(args.run_dir, "*.var")))
    if not vfiles:
        sys.exit(f"ERROR: no .var files in {args.run_dir}")

    # validate the target is a reference adenosine (checked against every well)
    for v in vfiles:
        c = counts_at(v, args.pos)
        if c is None:
            sys.exit(f"ERROR: position {args.pos} missing from {os.path.basename(v)}; "
                     f"cannot distinguish from zero editing.")
        if c["ref"] != "A":
            sys.exit(f"ERROR: position {args.pos} is '{c['ref']}', not a reference adenosine "
                     f"in {os.path.basename(v)}; check --pos.")
    print(f"# run={run}  ACTB target pos={args.pos} (validated ref=A)", file=sys.stderr)

    out = []
    for v in vfiles:
        well = os.path.splitext(os.path.basename(v))[0]
        c = counts_at(v, args.pos)
        A, G, ag = c["A"], c["G"], c["A"] + c["G"]
        edit = 100 * G / ag if ag else None
        other = c["C"] + c["T"]
        flags = []
        if c["reads_all"] < args.min_cov: flags.append("LOW_COV")
        if ag and other / (ag + other) > 0.02: flags.append("HIGH_NON_AG")
        if edit is None: flags.append("NO_A_OR_G")
        out.append([run, well, "ACTB", args.pos,
                    None if edit is None else round(edit, 4),
                    A, G, c["C"], c["T"], c["reads_all"], ";".join(flags)])

    w = csv.writer(sys.stdout)
    w.writerow(["run","well","site","pos","edit_pct","A","G","C","T","reads_all","flags"])
    w.writerows(out)

if __name__ == "__main__":
    main()

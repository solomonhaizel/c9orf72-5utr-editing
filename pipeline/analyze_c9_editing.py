#!/usr/bin/env python3
"""
A-to-I editing at the C9ORF72 5'UTR targets A74 and A88 from Genotyper .var output.

Editing = 100 * G / (A + G) at the target adenosine. This script reports raw
per-well editing and coverage with explicit validation; it does NOT subtract a
baseline (that is done per replicate against the matched same-run no-guide in
make_editing_table.py, so the run/rep pairing is explicit rather than inferred).

Validation performed (all failures raise or warn, never silently return 0):
  - ATGTGA start-stop motif must occur exactly once in the reference.
  - A88 ( = A74 + 14 ) must be a reference adenosine.
  - each target position must be present in every .var file.
  - coverage below --min-cov is reported and flagged, not dropped.
  - non-A/G reads at a target (C/T/indel share) are reported for QC.

Usage:
  python analyze_c9_editing.py RUN_DIR [--run NAME] [--min-cov 10000]
  RUN_DIR contains ref0.fa and the Plate1-*.var files.
"""
import argparse, csv, glob, os, sys

def load_ref(run_dir):
    seq = ""
    started = False
    with open(os.path.join(run_dir, "ref0.fa")) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if started:
                    break
                started = True
            elif started:
                seq += line
    return seq.upper()

def locate_sites(ref):
    hits = [i for i in range(len(ref)) if ref.startswith("ATGTGA", i)]
    if len(hits) != 1:
        sys.exit(f"ERROR: ATGTGA motif found {len(hits)} times (expected exactly 1).")
    a74 = hits[0] + 1
    a88 = a74 + 14
    if a88 > len(ref) or ref[a88-1] != "A":
        sys.exit(f"ERROR: A88 (pos {a88}) is not a reference adenosine "
                 f"(found '{ref[a88-1] if a88<=len(ref) else 'out-of-range'}').")
    return a74, a88

def counts_at(var_path, pos):
    for r in csv.DictReader(open(var_path), delimiter="\t"):
        if int(r["pos"]) == pos:
            return {b: int(r[b]) for b in "ACGT"} | {"reads_all": int(r["reads_all"])}
    return None   # position absent -> caller must handle explicitly

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--run", default=None, help="run/replicate label for the output")
    ap.add_argument("--min-cov", type=int, default=10000)
    args = ap.parse_args()
    run = args.run or os.path.basename(os.path.normpath(args.run_dir))

    ref = load_ref(args.run_dir)
    a74, a88 = locate_sites(ref)
    print(f"# run={run}  ref_len={len(ref)}  A74={a74}  A88={a88}", file=sys.stderr)

    out = []
    for var in sorted(glob.glob(os.path.join(args.run_dir, "*.var"))):
        well = os.path.splitext(os.path.basename(var))[0]
        for site, pos in (("A74", a74), ("A88", a88)):
            c = counts_at(var, pos)
            if c is None:
                sys.exit(f"ERROR: {site} (pos {pos}) missing from {well}.var; "
                         f"cannot distinguish from zero editing.")
            A, G, ag = c["A"], c["G"], c["A"] + c["G"]
            edit = 100 * G / ag if ag else None
            other = c["C"] + c["T"]
            flags = []
            if c["reads_all"] < args.min_cov: flags.append("LOW_COV")
            if ag and other / (ag + other) > 0.02: flags.append("HIGH_NON_AG")
            if edit is None: flags.append("NO_A_OR_G")
            out.append([run, well, site, pos,
                        None if edit is None else round(edit, 4),
                        A, G, c["C"], c["T"], c["reads_all"], ";".join(flags)])

    w = csv.writer(sys.stdout)
    w.writerow(["run","well","site","pos","edit_pct","A","G","C","T","reads_all","flags"])
    w.writerows(out)
    flagged = [r for r in out if r[-1]]
    if flagged:
        print(f"# {len(flagged)} well x site rows flagged (see 'flags' column)", file=sys.stderr)

if __name__ == "__main__":
    main()

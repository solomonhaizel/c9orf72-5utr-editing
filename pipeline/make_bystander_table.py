#!/usr/bin/env python3
"""
Build the amplicon-wide adenosine-scan table (bystander_by_rep.csv) from .var
output, with explicit condition and replicate IDs matching editing_by_rep.csv.

For every reference adenosine in the biological 1-98 window (16 positions,
coords ~30-96), records raw editing 100*G/(A+G) per condition per replicate.
Target sites A74 and A88 are flagged. This is a per-position raw scan (no
baseline subtraction); the No-guide and Scramble conditions are the controls.

'Amplicon-wide adenosine scan' is accurate: only the 16 reference adenosines in
the sequenced window are scored, not every base.
"""
import csv, os, sys
OFF = 33  # biological coord = ref position - 33
VAR = "data/raw_var"; CSV = "data/raw_perposition"

# condition -> list of (rep, run, source, guide-well/sample)
COND = {
 "A74 0.3 uM": [(1,"run1","csv","A743P"),(2,"run2","var","run2/Plate1-A2"),
               (3,"run3","var","run3/Plate1-A2"),(4,"run4","var","run4/Plate1-A2")],
 "A74 1 uM":   [(1,"run1","csv","A741uP"),(2,"run2","var","run2/Plate1-A1"),
               (3,"run3","var","run3/Plate1-A1"),(4,"run4","var","run4/Plate1-A1")],
 "A88 0.3 uM": [(1,"run1","csv","A883P"),(2,"run6","var","run6/Plate1-A2"),
               (3,"run5","var","run5/Plate1-A3"),(4,"run5","var","run5/Plate1-A4")],
 "A88 1 uM":   [(1,"run1","csv","A881uP"),(2,"run6","var","run6/Plate1-A1"),
               (3,"run5","var","run5/Plate1-A1"),(4,"run5","var","run5/Plate1-A2")],
 "No guide":  [(1,"run1","csv","UTP"),(2,"run2","var","run2/Plate1-A4"),
               (3,"run6","var","run6/Plate1-A7"),(4,"run5","var","run5/Plate1-A7")],
 "Scramble":  [(1,"run2","var","run2/Plate1-A3"),(2,"run3","var","run3/Plate1-A3"),
               (3,"run6","var","run6/Plate1-A6"),(4,"run5","var","run5/Plate1-A5")],
}

def a_positions():
    for line in open(f"{VAR}/run4/ref0.fa"):
        if not line.startswith(">"): 
            ref = line.strip().upper(); break
    return sorted({p-OFF for p in range(1,len(ref)+1) if ref[p-1]=="A" and 1<=p-OFF<=98})

def var_scan(path, apos):
    rows={int(r["pos"]): r for r in csv.DictReader(open(path), delimiter="\t")}
    out={}
    for bio in apos:
        r=rows.get(bio+OFF)
        if r:
            A,G=int(r["A"]),int(r["G"]); out[bio]=100*G/(A+G) if A+G else None
    return out

def csv_scan(sample, apos):
    d={(x["sample"],int(x["position"])):x for x in csv.DictReader(open(f"{CSV}/C9-Main.csv"))}
    out={}
    for bio in apos:
        r=d.get((sample,bio))
        if r:
            A,G=float(r["A_count"]),float(r["G_count"]); out[bio]=100*G/(A+G) if A+G else None
    return out

def main():
    apos=a_positions()
    out=[]
    for cond, reps in COND.items():
        for rep, run, src, key in reps:
            scan = var_scan(f"{VAR}/{key}.var", apos) if src=="var" else csv_scan(key, apos)
            for bio, e in scan.items():
                if e is not None:
                    out.append([cond, rep, run, bio, round(e,4), "target" if bio in (74,88) else "other"])
    with open("data/bystander_by_rep.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["condition","rep","run","position","edit_pct","site_type"]); w.writerows(out)
    print(f"wrote data/bystander_by_rep.csv ({len(out)} rows, {len(apos)} A positions)")

if __name__ == "__main__":
    main()

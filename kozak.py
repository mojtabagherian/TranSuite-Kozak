#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Kozak enrichment on CHANGED transcripts only.
- Filters dest-compare.csv to rows where PTC status changes between Reference and Revised.
- Extracts -3..+4 context around the Ref and Revised start codons.
- Scores strict Kozak match using RNNATGGV (R=A/G, V=A/C/G).
- Outputs a CSV with per-transcript results.
"""

import argparse, csv, re
from collections import defaultdict

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compare_csv", required=True, help="CSV with per-transcript PTC flags")
    ap.add_argument("--ref_gtf", required=True, help="Reference (Araport11) GTF")
    ap.add_argument("--rev_gtf", required=True, help="Revised (TranSuite) GTF")
    ap.add_argument("--tx_fasta", required=True, help="Transcriptome FASTA keyed by transcript_id")
    ap.add_argument("--tx_col", default="transcript_id", help="Column for transcript IDs in compare CSV")
    ap.add_argument("--ref_col", default="is_PTC50nt_ref", help="Column for Reference PTC flag")
    ap.add_argument("--revise_col", default="is_PTC50nt_revise", help="Column for Revised PTC flag")
    ap.add_argument("--out_prefix", default="kozak_changed", help="Output file prefix")
    return ap.parse_args()

def read_fasta(path):
    seqs = {}
    with open(path) as f:
        tid, buf = None, []
        for line in f:
            if line.startswith(">"):
                if tid: seqs[tid] = "".join(buf).upper()
                tid = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if tid: seqs[tid] = "".join(buf).upper()
    return seqs

def gtf_first_cds_start_on_transcript(gtf_path):
    exon_map, cds_map, strand = defaultdict(list), defaultdict(list), {}
    with open(gtf_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            chrom, src, ftype, start, end, score, st, phase, attrs = parts
            start, end = int(start), int(end)
            tid = None
            for a in attrs.split(";"):
                a = a.strip()
                if a.startswith("transcript_id"):
                    try: tid = a.split('"')[1]
                    except: tid = a.split()[-1]
            if not tid: continue
            strand[tid] = st
            if ftype == "exon": exon_map[tid].append((start,end))
            elif ftype == "CDS": cds_map[tid].append((start,end))
    tx_cds_starts = {}
    for tid in cds_map:
        if tid not in exon_map: continue
        exons = sorted(exon_map[tid], key=lambda x: x[0])
        cds   = sorted(cds_map[tid],  key=lambda x: x[0])
        st = strand.get(tid, "+")
        if st == "-": exons.reverse(); cds.reverse()
        tpos = 0; exon_blocks=[]
        for s,e in exons:
            length = e-s+1
            exon_blocks.append((s,e,tpos+1,tpos+length))
            tpos+=length
        g_start = cds[0][0] if st=="+" else cds[0][1]
        for gs,ge,ts,te in exon_blocks:
            if gs<=g_start<=ge:
                tx_start = ts+(g_start-gs) if st=="+" else ts+(ge-g_start)
                tx_cds_starts[tid]=tx_start
                break
    return tx_cds_starts

def parse_bool(v):
    if v is None: return None
    s=str(v).lower()
    if s in ("1","true","yes","y"): return True
    if s in ("0","false","no","n"): return False
    return None

def read_changed(compare_csv,tx_col,ref_col,revise_col):
    ids=[]
    with open(compare_csv,newline="") as f:
        rdr=csv.DictReader(f)
        for row in rdr:
            tid=row[tx_col].strip()
            p=parse_bool(row.get(ref_col))
            r=parse_bool(row.get(revise_col))
            if p!=r: ids.append(tid)
    return ids

def get_8mer(seq,start_1based):
    if not start_1based: return None
    i=start_1based-1
    if i-3<0 or i+5>len(seq): return None
    return seq[i-3:i+5]

def kozak_match(mer8):
    return bool(mer8 and re.fullmatch(r"[AG]..ATGG[ACG]", mer8))

def main():
    args=parse_args()
    changed_ids=read_changed(args.compare_csv,args.tx_col,args.ref_col,args.revise_col)
    ref_starts=gtf_first_cds_start_on_transcript(args.ref_gtf)
    rev_starts=gtf_first_cds_start_on_transcript(args.rev_gtf)
    seqs=read_fasta(args.tx_fasta)
    out_csv=args.out_prefix+"_results.csv"
    with open(out_csv,"w",newline="") as out:
        w=csv.writer(out)
        w.writerow(["transcript_id","ref_start","ref_8mer","ref_kozak",
                    "rev_start","rev_8mer","rev_kozak"])
        for tid in changed_ids:
            seq=seqs.get(tid)
            if not seq: continue
            rs,vs=ref_starts.get(tid),rev_starts.get(tid)
            ref8,rev8=get_8mer(seq,rs),get_8mer(seq,vs)
            rk,vk=kozak_match(ref8),kozak_match(rev8)
            w.writerow([tid,rs,ref8 or "NA",rk,vs,rev8 or "NA",vk])
    print(f"Done. Results in {out_csv}")

if __name__=="__main__": 
    main()
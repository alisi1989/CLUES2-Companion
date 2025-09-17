#!/usr/bin/env python3

#Example of usage for relate
#python CLUES2Companion_extract_freqs.py \
 # --engine Relate \
 # --chr 2 \
 # --vcf_prefix example/Finnish \
 # --out_prefix Finnish \
 # --start 135000000 --end 136000000 \
 # --haps output_CLUES2Companion-Relate/phase1/Finnish_chr2/Finnish_PIF_chr2.haps.gz \
 # --mut  output_CLUES2Companion-Relate/phase1/Finnish_chr2/Finnish_GS+COAL_chr2.mut

#Example of usare for Singer:
#python CLUES2Companion_extract_freqs.py \
 # --engine Singer \
 # --chr 2 \
 # --vcf_prefix example/Finnish \
 # --out_prefix Finnish \
 # --start 135000000 --end 136000000

import argparse
import os, sys, gzip, re
from cyvcf2 import VCF

# ---------- Helpers ----------------------------------------------------------
def info(msg): print(f"\033[1;32m→ {msg}\033[0m")
def warn(msg): print(f"\033[33m[WARN] {msg}\033[0m", file=sys.stderr)
def error(msg): sys.exit(f"\033[31m[ERROR] {msg}\033[0m")

def op(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def norm_chr(s: str) -> str:
    return s[3:] if str(s).lower().startswith("chr") else str(s)

# ---------- SNP extraction ---------------------------------------------------
def extract_snps(vcf_gz, chr, start, end, out_prefix, work_dir):
    out_file = os.path.join(work_dir, f"{out_prefix}_SNPs_chr{chr}_{start}_{end}.txt")
    vcf = VCF(vcf_gz)
    chr_norm = norm_chr(chr)
    n = 0
    with open(out_file, "w") as out:
        for v in vcf:
            if norm_chr(v.CHROM) != chr_norm: continue
            if v.POS < start or v.POS > end: continue
            if v.ID and v.ID != ".":
                out.write(f"{v.ID}\t{v.POS}\n")
                n += 1
    info(f"Extracted {n} SNPs → {out_file}")
    if n == 0: warn("No SNPs extracted in region")
    return out_file

# ---------- Frequenze & Derived alleles (Relate) ------------------------------
def calc_freq_relate(snp_file, haps_file, mut_file, work_dir, out_prefix, chr, start, end):
    snps = [(rs,pos) for rs,pos in (ln.split()[:2] for ln in open(snp_file))]
    pos2rs = {int(pos): rs for rs,pos in snps}

    # HAPS dictionary
    haps_by_rs, haps_by_pos = {}, {}
    with op(haps_file) as f:
        for ln in f:
            p = ln.split()
            if len(p) < 6: continue
            try: pos = int(p[2])
            except ValueError: continue
            rsid, a1, a2, g = p[1], p[3].upper(), p[4].upper(), [int(x) for x in p[5:]]
            haps_by_pos[pos] = (a1,a2,g)
            if rsid and rsid!=".": haps_by_rs[rsid] = (a1,a2,g,pos)

    # MUT dictionary
    rs2_ancder, pos2_ancder = {}, {}
    allele_pair_re = re.compile(r"^[ACGT]/[ACGT]$", re.IGNORECASE)
    with op(mut_file) as f:
        next(f) # skip header
        for ln in f:
            fs = ln.strip().split(";")
            if len(fs) >= 11 and fs[3] != ".":
                rs = fs[3]; ancder = fs[10]; flip = fs[7]
                if allele_pair_re.match(ancder):
                    anc, der = ancder.split("/")
                    if flip == "1": anc, der = der, anc
                    rs2_ancder[rs] = der.upper()
                    continue
            # fallback POS-based
            toks = ln.split()
            der=None; pos=None
            for t in toks:
                if allele_pair_re.match(t): _, der = t.split("/")
            for t in toks:
                if t.isdigit(): pos=int(t); break
            if der and pos: pos2_ancder[pos]=der

    # Output
    freq_file = os.path.join(work_dir, f"{out_prefix}_Frequency_chr{chr}_{start}_{end}.txt")
    with open(freq_file,"w") as out:
        out.write("rsID\tPOS\tallele_der\tfreq_der\n")
        written=0
        for rs,pos in snps:
            pos=int(pos)
            rec = haps_by_rs.get(rs) or haps_by_pos.get(pos)
            if not rec: continue
            if len(rec)==4: a1,a2,g,_=rec
            else: a1,a2,g=rec
            der = rs2_ancder.get(rs) or pos2_ancder.get(pos)
            if not der or der not in (a1,a2): continue
            bits = g if der==a2 else [1-b for b in g]
            freq = sum(bits)/len(bits)
            dfile = os.path.join(work_dir,f"{out_prefix}_Derived_{rs}.txt")
            with open(dfile,"w") as d: d.write("\n".join(map(str,bits)))
            out.write(f"{rs}\t{pos}\t{der}\t{freq:.4f}\n")
            written+=1
    info(f"✓ Frequency table → {freq_file} ({written} SNPs)")
    return freq_file

# ---------- Frequenze (Singer) ------------------------------------------------
def calc_freq_singer(snp_file, vcf_gz, work_dir, out_prefix, chr, start, end):
    freq_file = os.path.join(work_dir, f"{out_prefix}_Frequency_chr{chr}_{start}_{end}.txt")
    snps = {ln.split()[0] for ln in open(snp_file)}
    chr_norm = norm_chr(chr)
    vcf = VCF(vcf_gz)
    n=0
    with open(freq_file,"w") as out:
        out.write("rsID\tPOS\tallele_ALT\tfreq_ALT\n")
        for v in vcf:
            if not v.ID or v.ID=="." or v.ID not in snps: continue
            if norm_chr(v.CHROM)!=chr_norm: continue
            if v.POS<start or v.POS>end: continue
            alleles=[a for gt in v.genotypes for a in gt[:2] if a in (0,1)]
            if not alleles: continue
            alt_freq=sum(1 for a in alleles if a==1)/len(alleles)
            out.write(f"{v.ID}\t{v.POS}\t{v.ALT[0]}\t{alt_freq:.4f}\n")
            n+=1
    info(f"✓ Frequency table → {freq_file} ({n} SNPs)")
    return freq_file

# ---------- Main --------------------------------------------------------------
if __name__=="__main__":
    ap=argparse.ArgumentParser(description="CLUES2Companion SNP/Derived/Freq extractor")
    ap.add_argument("--engine",choices=["Relate","Singer"],required=True)
    ap.add_argument("--chr",required=True)
    ap.add_argument("--vcf_prefix",required=True,help="Prefix of VCF (without _chrN.vcf.gz)")
    ap.add_argument("--out_prefix",required=True)
    ap.add_argument("--start",type=int,required=True)
    ap.add_argument("--end",type=int,required=True)
    ap.add_argument("--haps",help="HAPS file (Relate)")
    ap.add_argument("--mut",help="MUT file (Relate)")
    args=ap.parse_args()

    work_base=os.path.join(os.getcwd(),f"output_CLUES2Companion-{args.engine}")
    run_id=f"{args.out_prefix}_chr{args.chr}"
    work_dir=os.path.join(work_base,"extract_only",run_id)
    os.makedirs(work_dir,exist_ok=True)

    vcf_gz=f"{args.vcf_prefix}_chr{args.chr}.vcf.gz"
    if not os.path.exists(vcf_gz): error(f"Missing {vcf_gz}")

    snp_file=extract_snps(vcf_gz,args.chr,args.start,args.end,args.out_prefix,work_dir)

    if args.engine=="Relate":
        if not args.haps or not args.mut: error("Need --haps and --mut for Relate mode")
        calc_freq_relate(snp_file,args.haps,args.mut,work_dir,args.out_prefix,args.chr,args.start,args.end)
    else:
        calc_freq_singer(snp_file,vcf_gz,work_dir,args.out_prefix,args.chr,args.start,args.end)

    info(f"✓ Extraction completed → {work_dir}")

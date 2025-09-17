#!/usr/bin/env python3
import argparse
import os
import sys
import subprocess
import gzip
from cyvcf2 import VCF

# ---------- Helpers ----------------------------------------------------------

def info(msg): print(f"\033[1;32m→ {msg}\033[0m")
def warn(msg): print(f"\033[33m[WARN] {msg}\033[0m", file=sys.stderr)
def error(msg): sys.exit(f"\033[31m[ERROR] {msg}\033[0m")

def abs_path(p: str) -> str:
    if not p: return ""
    return os.path.abspath(os.path.expanduser(p))

def run_cmd(cmd, cwd=None):
    info(f"Running: {' '.join(map(str,cmd))}")
    try:
        subprocess.run(cmd, check=True, cwd=cwd)
    except subprocess.CalledProcessError:
        error(f"Command failed: {' '.join(map(str,cmd))}")

import random

def run_phase1_relate(args):
    out_prefix = args.out_prefix
    chr        = args.chr
    start_bp   = args.start
    end_bp     = args.end

    work_base = os.path.join(os.getcwd(), "output_CLUES2Companion-Relate")
    run_id    = f"{out_prefix}_chr{chr}"
    work_dir  = os.path.join(work_base, "phase1", run_id)
    os.makedirs(work_dir, exist_ok=True)
    info(f"Working dir: {work_dir}")

    relate_root = os.path.join(os.getcwd(), "Relate-Linux")
    relate_bin  = os.path.join(relate_root, "bin", "Relate")
    if not os.path.exists(relate_bin):
        error(f"Relate binary not found at {relate_bin}")

    # input files
    vcf_prefix = f"{args.vcf_prefix}_chr{chr}"
    haps_path  = os.path.join(work_dir, f"{out_prefix}_chr{chr}.haps")
    sample_path= os.path.join(work_dir, f"{out_prefix}_chr{chr}.sample")

    # internal prefixes
    pif_base   = os.path.join(work_dir, f"{out_prefix}_PIF_chr{chr}")
    gs_base    = os.path.join(work_dir, f"{out_prefix}_GS_chr{chr}")
    gs_coal    = os.path.join(work_dir, f"{out_prefix}_GS+COAL_chr{chr}")
    coal_base  = os.path.join(work_dir, f"{out_prefix}_EPS4COAL_chr{chr}")

    # 4. Convert VCF → haps/sample
    run_cmd([os.path.join(relate_root,"bin/RelateFileFormats"),
             "--mode","ConvertFromVcf",
             "-i", vcf_prefix,
             "--haps", haps_path,
             "--sample", sample_path,
             "-o", os.path.join(work_dir,f"{out_prefix}_chr{chr}")])

    # 5. PrepareInputFiles
    rand_seed = random.randint(1, 1_000_000_000)
    run_cmd([os.path.join(relate_root,"scripts/PrepareInputFiles/PrepareInputFiles.sh"),
             "--haps", haps_path,
             "--sample", sample_path,
             "--ancestor", args.ancestor,
             "--mask", args.mask,
             "--poplabels", args.poplabels,
             "--seed", str(rand_seed),
             "-o", pif_base])

    pif_haps   = f"{pif_base}.haps.gz"
    pif_sample = f"{pif_base}.sample.gz"
    pif_annot  = f"{pif_base}.annot"

    # 6. Relate run-1
    rand_seed = random.randint(1, 1_000_000_000)
    run_cmd([
        relate_bin,
        "--mode","All",
        "--haps", pif_haps,
        "--sample", pif_sample,
        "--map", args.map,
        "--annot", pif_annot,
        "-N", str(args.N),
        "-m", str(args.m),
        "-o", f"{out_prefix}_GS_chr{chr}_run1",
        "--seed", str(rand_seed)
    ], cwd=work_dir)

    # 7. EstimatePopulationSize
    rand_seed = random.randint(1, 1_000_000_000)
    run_cmd([os.path.join(relate_root,"scripts/EstimatePopulationSize/EstimatePopulationSize.sh"),
             "-i", f"{out_prefix}_GS_chr{chr}_run1",
             "-o", f"{out_prefix}_EPS4COAL_chr{chr}",
             "--noanc", "1",
             "--poplabels", args.poplabels,
             "-m", str(args.m),
             "--years_per_gen","28",
             "--seed", str(rand_seed)
    ], cwd=work_dir)

    # 8. Re-estimate branch lengths
    rand_seed = random.randint(1, 1_000_000_000)
    run_cmd([os.path.join(relate_root,"scripts/SampleBranchLengths/ReEstimateBranchLengths.sh"),
             "-i", f"{out_prefix}_GS_chr{chr}_run1",
             "-o", f"{out_prefix}_GS+COAL_chr{chr}",
             "-m", str(args.m),
             "--coal", f"{out_prefix}_EPS4COAL_chr{chr}.coal",
             "--seed", str(rand_seed)
    ], cwd=work_dir)

    info("✓ Relate pipeline completed")


# ---------- Phase-1 Singer ---------------------------------------------------

def run_phase1_singer(args):
    out_prefix = args.out_prefix
    chr        = args.chr
    start_bp   = args.start
    end_bp     = args.end

    work_base = os.path.join(os.getcwd(), "output_CLUES2Companion-Singer")
    run_id    = f"{out_prefix}_chr{chr}"
    work_dir  = os.path.join(work_base, "phase1", run_id)
    os.makedirs(work_dir, exist_ok=True)
    info(f"Working dir: {work_dir}")

    singer_exec = os.path.join(os.getcwd(), "Singer-Linux", "singer_master")
    if not os.path.exists(singer_exec):
        error(f"SINGER executable not found: {singer_exec}")

    vcf_plain = f"{args.vcf_prefix}_chr{chr}"
    vcf_gz    = f"{args.vcf_prefix}_chr{chr}.vcf.gz"

    cmd = [singer_exec, "-vcf", vcf_plain, "-output", os.path.join(work_dir, out_prefix),
           "-start", str(start_bp), "-end", str(end_bp)]
    if args.Ne:      cmd += ["-Ne", str(args.Ne)]
    if args.ratio:   cmd += ["-ratio", str(args.ratio)]
    if args.recomb:  cmd += ["-recomb_map", args.recomb]
    if args.mutmap:  cmd += ["-mut_map", args.mutmap]
    if args.n:       cmd += ["-n", str(args.n)]
    if args.thin:    cmd += ["-thin", str(args.thin)]
    if args.polar:   cmd += ["-polar", str(args.polar)]
    if args.seed:    cmd += ["-seed", str(args.seed)]
    if args.mu:    cmd += ["-m", str(args.mu)]

    run_cmd(cmd)

    # Convert to tskit
    convert_exec = os.path.join(os.getcwd(), "Singer-Linux", "convert_to_tskit")
    if os.path.exists(convert_exec):
        run_cmd([convert_exec, "-input", os.path.join(work_dir,out_prefix),
                 "-output", os.path.join(work_dir,out_prefix),
                 "-start","0","-end", str(args.n or 100)])
    else:
        warn("convert_to_tskit not found, skipping")

    info("✓ Singer pipeline completed")
    
# ---------- Frequenze & Derived alleles (versione consolidata) ---------------

def calc_freq_relate(snp_file, haps_file, mut_file, work_dir, out_prefix, chr, start, end):
    """
    Calcola DER e frequenze per Relate:
    - preferisce la join per rsID (come nello script Bash),
    - fallback a join per POS,
    - supporta .mut o .mut.gz e .haps.gz.
    """
    import re, gzip

    def op(path):
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

    # --- SNP list (rs, pos) --------------------------------------------------
    snps = []
    pos2rs = {}
    with open(snp_file) as f:
        for ln in f:
            rs, pos = ln.split()[:2]
            pos = int(pos)
            snps.append((rs, pos))
            pos2rs[pos] = rs

    # --- HAPS index: by rsID e by POS ---------------------------------------
    # Formato .haps(.gz): CHR ID POS A1 A2 H1 H2 ...
    haps_by_rs = {}
    haps_by_pos = {}
    with op(haps_file) as f:
        for ln in f:
            p = ln.split()
            if len(p) < 6:
                continue
            try:
                pos = int(p[2])
            except ValueError:
                continue
            rsid = p[1]
            a1, a2 = p[3].upper(), p[4].upper()
            try:
                g = [int(x) for x in p[5:]]  # 0→A1, 1→A2
            except ValueError:
                continue
            haps_by_pos[pos] = (a1, a2, g)
            if rsid and rsid != ".":
                haps_by_rs[rsid] = (a1, a2, g, pos)

    # --- MUT index: prova parsing "classico" (a;...;rs;...;flip;...;ANC/DER) -
    #               + parser resiliente come fallback
    rs2_anc_der = {}   # rs -> (ANC, DER, flip_used)
    pos2_anc_der = {}  # pos -> (ANC, DER)

    allele_pair_re = re.compile(r"^[ACGT]/[ACGT]$", re.IGNORECASE)

    with op(mut_file) as f:
        header_seen = False
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            if not header_seen:
                header_seen = True
                # lo script Bash fa next(f), quindi qui saltiamo solo la prima riga
                # senza ispezionarla
                if s.startswith("#") or s.startswith("index"):
                    continue
            # 1) parser "classico" a ';'
            parts = s.split(";")
            if len(parts) >= 11:
                rs = parts[3]
                anc_der = parts[10].strip().upper() if parts[10] else ""
                flip = parts[7].strip() if len(parts) > 7 else "0"
                if rs and rs != "." and allele_pair_re.match(anc_der):
                    anc, der = anc_der.split("/")
                    if flip == "1":
                        anc, der = der, anc
                    rs2_anc_der[rs] = (anc, der, True)
                    # prova anche a riempire pos come fallback (primo int plausibile)
                    for t in parts:
                        if t.isdigit():
                            pos = int(t)
                            pos2_anc_der[pos] = (anc, der)
                            break
                    continue  # passa alla prossima riga .mut

            # 2) fallback: token whitespace e ricerca 'A/C' + primo intero come POS
            toks = s.replace("\t", " ").split()
            anc, der, pos = None, None, None
            for t in toks:
                if allele_pair_re.match(t.upper()):
                    anc, der = t.upper().split("/")
                    break
            for t in toks:
                if t.isdigit():
                    pos = int(t)
                    break
            if anc and der and pos is not None:
                pos2_anc_der[pos] = (anc, der)

    # --- Output --------------------------------------------------------------
    out_freq = os.path.join(work_dir, f"{out_prefix}_Frequency_chr{chr}_{start}_{end}.txt")

    miss_haps = miss_mut = mismatch_alleles = 0
    written = total = 0

    with open(out_freq, "w") as out:
        out.write("rsID\tPOS\tallele_der\tfreq_der\n")
        for rs, pos in snps:
            total += 1

            # 1) recupero alleli e genotipi dall'HAPS
            rec = haps_by_rs.get(rs) or haps_by_pos.get(pos)
            if rec is None:
                miss_haps += 1
                continue
            if len(rec) == 4:
                a1, a2, g, pos_haps = rec
            else:
                a1, a2, g = rec
                pos_haps = pos

            # 2) recupero ANC/DER dal MUT (preferenza rs, poi pos)
            tpl = rs2_anc_der.get(rs)
            if tpl:
                anc, der, _ = tpl
            else:
                tpl2 = pos2_anc_der.get(pos)
                if not tpl2:
                    miss_mut += 1
                    continue
                anc, der = tpl2

            der = der.upper()
            a1, a2 = a1.upper(), a2.upper()

            if der not in (a1, a2):
                # a volte REF/ALT nell'HAPS differiscono (liftover/strand). Senza ANC certo, meglio
                # non forzare: segnalo mismatch e salto.
                mismatch_alleles += 1
                continue

            # 3) polarizzo: 1 == DER
            bits = g if der == a2 else [1 - b for b in g]
            freq = sum(bits) / len(bits) if bits else 0.0

            rs_out = rs if rs and rs != "." else f"pos{pos_haps}"
            dpath = os.path.join(work_dir, f"{out_prefix}_Derived_{rs_out}.txt")
            with open(dpath, "w") as d:
                d.write("\n".join(map(str, bits)))

            out.write(f"{rs_out}\t{pos_haps}\t{der}\t{freq:.4f}\n")
            written += 1

    info(f"✓ Frequency table → {out_freq}")
    info(f"  Stats: target={total}, written={written}, no_haps={miss_haps}, no_mut={miss_mut}, allele_mismatch={mismatch_alleles}")
    if written == 0:
        warn("Tabella vuota: probabile disallineamento rsID tra .mut e .haps oppure alleli non corrispondenti (liftover/strand).")
    return out_freq




def calc_freq_singer(snp_file, vcf_gz, work_dir, out_prefix, chr, start, end):
    freq_file = os.path.join(work_dir, f"{out_prefix}_Frequency_chr{chr}_{start}_{end}.txt")
    vcf = VCF(vcf_gz)

    snps = {ln.split()[0] for ln in open(snp_file)}
    def norm_chr(s): return s[3:] if str(s).lower().startswith("chr") else str(s)
    chr_norm = norm_chr(chr)

    with open(freq_file,"w") as out:
        out.write("rsID\tPOS\tallele_ALT\tfreq_ALT\n")
        n_written = 0
        for v in vcf:
            if not v.ID or v.ID=="." or v.ID not in snps: continue
            if norm_chr(v.CHROM) != chr_norm: continue
            if v.POS < start or v.POS > end: continue
            alleles = [a for pair in v.genotypes for a in pair[:2] if a in (0,1)]
            if not alleles: continue
            alt_freq = sum(1 for a in alleles if a==1) / len(alleles)
            alt_allele = v.ALT[0] if v.ALT else "."
            out.write(f"{v.ID}\t{v.POS}\t{alt_allele}\t{alt_freq:.4f}\n")
            n_written += 1

    info(f"✓ Frequency table → {freq_file}")
    return freq_file

# ---------- SNP extraction (common) ------------------------------------------

def extract_snps(args, work_dir, chr, start_bp, end_bp, out_prefix):
    vcf_gz = f"{args.vcf_prefix}_chr{chr}.vcf.gz"
    snp_file = os.path.join(work_dir, f"{out_prefix}_SNPs_chr{chr}_{start_bp}_{end_bp}.txt")

    vcf = VCF(vcf_gz)
    count = 0
    with open(snp_file,"w") as out:
        for v in vcf:
            if v.CHROM.replace("chr","") != str(chr): continue
            if v.POS < start_bp or v.POS > end_bp: continue
            if v.ID and v.ID!=".":
                out.write(f"{v.ID}\t{v.POS}\n")
                count += 1
    info(f"Extracted {count} SNPs → {snp_file}")
    return snp_file

# ---------- Main --------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CLUES2Companion Phase-1 CLI")
    parser.add_argument("--engine", choices=["Relate","Singer"], required=True, help="Inference engine")
    parser.add_argument("--chr", required=True, help="Chromosome (e.g. 2, X) [Relate/Singer]")
    parser.add_argument("--vcf_prefix", required=True, help="Prefix of VCF (without _chrN.vcf/.vcf.gz)[Relate/Singer]")
    parser.add_argument("--out_prefix", required=True, help="Output prefix name [Relate/Singer]")
    parser.add_argument("--start", type=int, required=True, help="Start bp of region [Relate/Singer]")
    parser.add_argument("--end", type=int, required=True, help="End bp of region [Relate/Singer]")
    # Relate params
    parser.add_argument("--poplabels", required=False, help="Population labels file (*.poplabels) [Relate]")
    parser.add_argument("-m", required=False, default="1.25e-8", help="Mutation rate (Relate)")
    parser.add_argument("-N", required=False, default="30000", help="Ne (Relate)")
    # Singer extras
    parser.add_argument("--mu", type=float, help="Mutation rate (Singer)")
    parser.add_argument("--Ne", type=int, help="Effective population size (Singer)")
    parser.add_argument("--ratio", type=float, help="Recomb/mutation ratio (Singer)")
    parser.add_argument("--recomb", help="Recombination map (Singer)")
    parser.add_argument("--mutmap", help="Mutation map (Singer)")
    parser.add_argument("-n", type=int, help="MCMC samples (Singer)")
    parser.add_argument("--thin", type=int, help="Thinning interval (Singer)")
    parser.add_argument("--polar", type=float, help="Site flip probability (Singer)")
    parser.add_argument("--seed", type=int, help="Random seed (Singer)")
    args = parser.parse_args()



if args.engine == "Relate":
    required_dir = os.path.join(os.getcwd(), "required_files")
    
    args.ancestor = os.path.join(required_dir, "ancestor", f"homo_sapiens_ancestor_chr{args.chr}.fa")
    args.mask     = os.path.join(required_dir, "mask", f"PilotMask_chr{args.chr}.fasta")
    args.map      = os.path.join(required_dir, "map", f"genetic_map_chr{args.chr}.txt")
    args.poplabels = abs_path(args.poplabels)

    for f in [args.ancestor, args.mask, args.map, args.poplabels]:
        if not os.path.exists(f):
            error(f"Missing required file: {f}")

    run_phase1_relate(args)   # <-- qui deve chiamare relate
else:  # Singer
    run_phase1_singer(args)


# --- step comuni: SNP extraction e frequenze ---
work_base = os.path.join(os.getcwd(), f"output_CLUES2Companion-{args.engine}")
run_id    = f"{args.out_prefix}_chr{args.chr}"
work_dir  = os.path.join(work_base,"phase1",run_id)

snp_file = extract_snps(args, work_dir, args.chr, args.start, args.end, args.out_prefix)

if args.engine == "Relate":
    mut_file  = os.path.join(work_dir, f"{args.out_prefix}_GS+COAL_chr{args.chr}.mut")
    if not os.path.exists(mut_file):
        if os.path.exists(mut_file + ".gz"):
            mut_file = mut_file + ".gz"
        else:
            error("Missing mut file for Relate (.mut o .mut.gz)")

    haps_file = os.path.join(work_dir, f"{args.out_prefix}_PIF_chr{args.chr}.haps.gz")
    if not os.path.exists(haps_file):
        # alcune pipeline salvano anche non compresso
        alt = haps_file[:-3]
        if os.path.exists(alt):
            haps_file = alt
        else:
            error("Missing haps file for Relate (.haps o .haps.gz)")

    calc_freq_relate(snp_file, haps_file, mut_file, work_dir,
                     args.out_prefix, args.chr, args.start, args.end)
else:
    vcf_gz = f"{args.vcf_prefix}_chr{args.chr}.vcf.gz"
    if not os.path.exists(vcf_gz): error(f"Missing VCF: {vcf_gz}")
    calc_freq_singer(snp_file, vcf_gz, work_dir,
                     args.out_prefix, args.chr, args.start, args.end)

info(f"✓ Phase-1 completed — outputs in {work_dir}")







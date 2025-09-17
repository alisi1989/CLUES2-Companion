#!/usr/bin/env python3
import argparse, os, sys, subprocess, random, glob, gzip, shlex

# ============== Helpers ======================================================

def info(msg): print(f"\033[1;32m→ {msg}\033[0m")
def warn(msg): print(f"\033[33m[WARN] {msg}\033[0m", file=sys.stderr)
def error(msg): sys.exit(f"\033[31m[ERROR] {msg}\033[0m")

def run_cmd(cmd, cwd=None):
    if isinstance(cmd, (list, tuple)):
        cmd_list = list(map(str, cmd))
    else:
        # allow string commands too
        cmd_list = shlex.split(cmd)
    info(f"Running: {' '.join(cmd_list)}")
    try:
        subprocess.run(cmd_list, check=True, cwd=cwd)
    except subprocess.CalledProcessError:
        error(f"Command failed: {' '.join(cmd_list)}")

def abs_path(p: str) -> str:
    if not p: return ""
    return os.path.abspath(os.path.expanduser(p))

# ============== Core =========================================================

def auto_phase1_dir(base_dir, engine, out_prefix, chr_):
    """
    Costruisce la cartella di Phase-1 coerente con Phase-1 CLI:
    output_CLUES2Companion-Relate/phase1/<OUT>_chrN
    output_CLUES2Companion-Singer/phase1/<OUT>_chrN
    """
    base_name = f"output_CLUES2Companion-{engine}"
    run_id    = f"{out_prefix}_chr{chr_}"
    return os.path.join(base_dir, base_name, "phase1", run_id)

def ensure_file(path, msg=None):
    if not os.path.isfile(path):
        error(msg or f"Missing file: {path}")

def ensure_dir(path, msg=None):
    if not os.path.isdir(path):
        error(msg or f"Missing directory: {path}")

def parse_timebins(tb_str):
    if not tb_str: return []
    # accetta "200 300" o "200,300"
    parts = [x for x in tb_str.replace(",", " ").split() if x.strip()]
    return parts

# ============== Phase-2 Runner ==============================================

def run_phase2(args):
    base_dir  = os.getcwd()
    engine    = args.engine
    chr_      = args.chr
    out_pref  = args.out_prefix

    # WORK_BASE dipende dal ramo
    work_base = os.path.join(base_dir, f"output_CLUES2Companion-{engine}")
    work_dir  = os.path.join(work_base, "phase2", f"{out_pref}_chr{chr_}")
    os.makedirs(work_dir, exist_ok=True)

    print("\n          ╔═════════════════════════════════════════════════════════╗")
    print("          ║       🧬  PHASE 2 – (requires Phase-1 outputs)  🧬      ║")
    print("          ║   Please read the manual carefully before proceeding    ║")
    print("          ╚═════════════════════════════════════════════════════════╝\n")

    # --- Phase-1 directory (auto-detect + override)
    phase1_dir = args.phase1_dir or auto_phase1_dir(base_dir, engine, out_pref, chr_)
    phase1_dir = abs_path(phase1_dir)
    ensure_dir(phase1_dir, f"Phase-1 directory not found: {phase1_dir}")

    # --- freq & snps (comuni)
    freq_glob = os.path.join(phase1_dir, f"{out_pref}_Frequency_chr{chr_}_*.txt")
    snps_glob = os.path.join(phase1_dir, f"{out_pref}_SNPs_chr{chr_}_*.txt")

    freq_files = sorted(glob.glob(freq_glob))
    snp_files  = sorted(glob.glob(snps_glob))
    if not freq_files: error(f"Frequency file not found in {phase1_dir}")
    if not snp_files:  error(f"SNP coordinate file not found in {phase1_dir}")

    freq_file = freq_files[0]
    snp_file  = snp_files[0]
    info(f"Using frequency file: {freq_file}")
    info(f"Using SNPs file: {snp_file}\n")

    # --- tool paths
    relate_root = os.path.join(base_dir, "Relate-Linux")     # coerente con Phase-1 CLI
    clues_dir   = os.path.join(base_dir, "CLUES2")
    rtc_py      = os.path.join(clues_dir, "RelateToCLUES.py")
    inf_py      = os.path.join(clues_dir, "inference.py")
    stc_py      = os.path.join(clues_dir, "SingerToCLUES.py")

    ensure_file(inf_py, "CLUES2/inference.py not found")
    if engine == "Relate":
        ensure_file(rtc_py, "CLUES2/RelateToCLUES.py not found")
        if not os.path.exists(os.path.join(relate_root, "bin", "Relate")):
            error("Relate binary not found in ./Relate/bin")
    else:
        ensure_file(stc_py, "CLUES2/SingerToCLUES.py not found")

    # --- CLUES2 params
    time_bins_list = parse_timebins(args.timeBins)
    # output dirs
    tree_dir  = os.path.join(work_dir, f"{out_pref}_trees_chr{chr_}")
    times_dir = os.path.join(work_dir, f"{out_pref}_times_chr{chr_}")
    infer_dir = os.path.join(work_dir, f"{out_pref}_inference_chr{chr_}")
    os.makedirs(tree_dir, exist_ok=True)
    os.makedirs(times_dir, exist_ok=True)
    os.makedirs(infer_dir, exist_ok=True)

    # -------------------- RELATE BRANCH --------------------------------------
    if engine == "Relate":
        # Artifacts di Phase-1
        gs_prefix  = os.path.join(phase1_dir, f"{out_pref}_GS+COAL_chr{chr_}")
        anc_ok = os.path.isfile(gs_prefix + ".anc") or os.path.isfile(gs_prefix + ".anc.gz")
        if not anc_ok:
            error(f"*.anc file not found for prefix {gs_prefix}")

        coal_file = os.path.join(phase1_dir, f"{out_pref}_EPS4COAL_chr{chr_}.coal")
        ensure_file(coal_file, f"*.coal file not found: {coal_file}")

        # Derived files
        derived_glob = os.path.join(phase1_dir, f"{out_pref}_Derived_*.txt")
        derived_list = sorted(glob.glob(derived_glob))
        if not derived_list:
            error(f"Derived files not found in {phase1_dir}")

        # Step A – SampleBranchLengths
        if args.num_samples is None or args.num_samples <= 0:
            error("For Relate branch, --num_samples must be a positive integer")

        info("Step A: SampleBranchLengths")
        with open(snp_file) as f:
            for ln in f:
                ln = ln.strip()
                if not ln or ln.startswith("#") or ln.lower().startswith("rsid"):
                    continue
                rsid, pos = ln.split()[:2]
                seed = random.randint(1, 1_000_000_000)
                cmd = [
                    os.path.join(relate_root, "scripts", "SampleBranchLengths", "SampleBranchLengths.sh"),
                    "--input", gs_prefix,
                    "--output", os.path.join(tree_dir, rsid),
                    "--first_bp", pos, "--last_bp", pos,
                    "--format", "n",
                    "--num_samples", str(args.num_samples),
                    "--coal", coal_file,
                    "-m", str(args.mu),
                    "--seed", str(seed)
                ]
                try:
                    run_cmd(cmd)
                except SystemExit:
                    warn(f"BranchLengths failed for {rsid}")
        print("✓ Step A completed")

        # Step B – RelateToCLUES
        info("Step B: RelateToCLUES")
        for dfile in derived_list:
            rsid = os.path.basename(dfile).removeprefix(f"{out_pref}_Derived_").removesuffix(".txt")
            nwx  = os.path.join(tree_dir, f"{rsid}.newick")
            if os.path.isfile(nwx):
                cmd = [sys.executable, rtc_py,
                       "--DerivedFile", dfile,
                       "--RelateSamples", nwx,
                       "--out", os.path.join(times_dir, rsid)]
                try:
                    run_cmd(cmd)
                except SystemExit:
                    warn(f"RTC failed for {rsid}")
            else:
                warn(f"newick not found for {rsid}")
        print("✓ Step B completed")

        singer_ne = None  # not used in this branch

    # -------------------- SINGER BRANCH --------------------------------------
    else:
        # Intervallo e Ne
        if args.region_start is None or args.region_end is None:
            error("For Singer branch, you must provide --region_start and --region_end")
        if args.Ne is None or args.Ne <= 0:
            error("For Singer branch, you must provide --Ne (effective population size)")

        singer_start = int(args.region_start)
        singer_end   = int(args.region_end)
        singer_ne    = int(args.Ne)

        trees_path = phase1_dir
        ensure_dir(trees_path, f"TSKIT trees directory not found: {trees_path}")

        # Estrai gli SNP nell'intervallo
        snp_list = []
        with open(snp_file) as f:
            for ln in f:
                ln = ln.strip()
                if not ln or ln.startswith("#") or ln.lower().startswith("rsid"):
                    continue
                rs, pos = ln.split()[:2]
                pos = int(pos)
                if singer_start <= pos <= singer_end:
                    snp_list.append((rs, pos))

        if not snp_list:
            error(f"No SNPs found in interval {singer_start}-{singer_end}")

        info(f"Found {len(snp_list)} SNPs in interval {singer_start}-{singer_end}")

        # Freq file "ridotto" con solo questi SNP
        freq_one = os.path.join(work_dir, "__freq_selected.txt")
        with open(freq_file) as fin, open(freq_one, "w") as fout:
            header = fin.readline().rstrip("\n")
            fout.write(header + "\n")
            rs_set = {rs for rs,_ in snp_list}
            for ln in fin:
                if not ln.strip(): continue
                parts = ln.split()
                if parts and parts[0] in rs_set:
                    fout.write(ln)
        freq_file = freq_one  # d'ora in poi usiamo solo questo

        # Loop SingerToCLUES per ogni posizione
        for rs, pos in snp_list:
            info(f"SingerToCLUES for {rs} at POS={pos}")
            cmd = [sys.executable, stc_py,
                   "--position", str(pos),
                   "--tree_path", trees_path,
                   "--output", os.path.join(times_dir, rs)]
            try:
                run_cmd(cmd)
            except SystemExit:
                warn(f"SingerToCLUES failed for {rs} ({pos})")

    # -------------------- Step C – inference.py (comune) ---------------------
    info("Step C: inference.py")
    with open(freq_file) as f:
        header = f.readline()  # skip
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            if len(parts) < 4:
                warn(f"Skipping malformed line in freq file: {ln}")
                continue
            rs, pos, allele, freq = parts[0], parts[1], parts[2], parts[3]
            tfile = os.path.join(times_dir, f"{rs}_times.txt")
            if not os.path.isfile(tfile):
                warn(f"times not found for {rs}")
                continue

            cmd = [sys.executable, inf_py,
                   "--times", tfile,
                   "--popFreq", freq,
                   "--out", os.path.join(infer_dir, rs),
                   "--tCutoff", str(args.tcutoff),
                   "--df", str(args.df)]
            # Relate: usa coal; Singer: usa Ne
            if engine == "Relate":
                coal_file = os.path.join(phase1_dir, f"{out_pref}_EPS4COAL_chr{chr_}.coal")
                if os.path.isfile(coal_file):
                    cmd += ["--coal", coal_file]
            else:
                cmd += ["--N", str(singer_ne)]

            info(f"Using dominance coefficient h={args.h}")


            if args.CI is not None:
                cmd += ["--CI", str(args.CI)]
            cmd += ["--h", str(args.h)]
            if args.anc_samps and os.path.isfile(args.anc_samps):
                cmd += ["--ancientSamps", args.anc_samps]
            if args.anc_haps and os.path.isfile(args.anc_haps):
                cmd += ["--ancientHaps", args.anc_haps]
            if args.noTraj:
                cmd += ["--noAlleleTraj"]
            if time_bins_list:
                cmd += ["--timeBins"] + time_bins_list

            try:
                run_cmd(cmd)
            except SystemExit:
                warn(f"inference failed for {rs}")
    print("✓ Step C completed")

    # -------------------- Step D – merge -------------------------------------
    info("Step D: merge inference outputs")
    merged = os.path.join(work_dir, f"{out_pref}_merged_inference_chr{chr_}.tsv")

    # trova un primo *_inference.txt per header/numero epoche
    inf_list = sorted(glob.glob(os.path.join(infer_dir, "*_inference.txt")))
    if not inf_list:
        error(f"No inference files found in {infer_dir}")

    first_inf = inf_list[0]
    with open(first_inf) as f:
        hdr_line = f.readline().rstrip("\n")
    n_sel = sum(1 for c in hdr_line.split("\t") if c.startswith("SelectionMLE"))

    with open(merged, "w") as out:
        out.write("rsID\tPOS\tder_freq\tlogLR\t-log10(p)")
        for k in range(1, n_sel+1):
            out.write(f"\tSelectionMLE{k}\tCI{k}_lower\tCI{k}_upper")
        out.write("\n")

        # ri-apriamo il freq_file (ridotto per Singer; completo per Relate)
        with open(freq_file) as fin:
            _ = fin.readline()  # header
            for ln in fin:
                ln = ln.strip()
                if not ln:
                    continue
                parts = ln.split("\t")
                if len(parts) < 4:
                    parts = ln.split()
                if len(parts) < 4:
                    warn(f"Skipping malformed freq row: {ln}")
                    continue
                rs, pos, _allele, freq = parts[:4]
                inf_file = os.path.join(infer_dir, f"{rs}_inference.txt")
                if not os.path.isfile(inf_file):
                    warn(f"inference missing for {rs}")
                    continue

                # header & data
                with open(inf_file) as fi:
                    hdr = fi.readline().rstrip("\n").split("\t")
                    data = fi.readline().rstrip("\n").split("\t")

                def col_idx(needle):
                    for i, h in enumerate(hdr):
                        if needle in h:
                            return i
                    return -1

                loglr_idx = col_idx("logLR")
                logp_idx  = col_idx("-log10(")
                if loglr_idx < 0 or logp_idx < 0:
                    warn(f"Header mismatch in {inf_file}")
                    continue

                row = [rs, pos, freq, data[loglr_idx], data[logp_idx]]

                # CI opzionale
                ci_file = os.path.join(infer_dir, f"{rs}_CI.txt")
                ci_lo, ci_up = [], []
                if os.path.isfile(ci_file):
                    with open(ci_file) as cf:
                        next(cf, None)  # skip header
                        for cl in cf:
                            epoch, lo, up = cl.rstrip("\n").split("\t")[:3]
                            ci_lo.append(lo); ci_up.append(up)

                for k in range(1, n_sel+1):
                    sel_idx = col_idx(f"SelectionMLE{k}")
                    sel_val = data[sel_idx] if sel_idx >= 0 and sel_idx < len(data) else "NA"
                    lo = ci_lo[k-1] if ci_lo and len(ci_lo) >= k else "NA"
                    up = ci_up[k-1] if ci_up and len(ci_up) >= k else "NA"
                    row += [sel_val, lo, up]

                out.write("\t".join(row) + "\n")

    info(f"Merged file → {merged}")

    # -------------------- Plot integrato -------------------------------------
    info("Generating integrated plot")
    try:
        import pandas as pd, numpy as np, matplotlib.pyplot as plt
        from adjustText import adjust_text
        from matplotlib.ticker import FuncFormatter
        from matplotlib.transforms import offset_copy

        df = pd.read_csv(merged, sep="\t")
        if "-log10(p)" not in df.columns:
            warn("Plot skipped: '-log10(p)' column not found in merged TSV.")
            return

        df["p_raw"] = 10 ** (-df["-log10(p)"])
        def stars(p):
            if p < 1e-3: return "***"
            if p < 1e-2: return "**"
            if p < 5e-2: return "*"
            return ""

        epochs = [c for c in df.columns if c.startswith("SelectionMLE")]
        n_ep   = len(epochs)
        prefix = os.path.join(work_dir, f"{out_pref}_clues")

        for k in range(1, n_ep+1):
            sel  = df[f"SelectionMLE{k}"]
            loCI = df[f"CI{k}_lower"]
            upCI = df[f"CI{k}_upper"]

            # yerr come nel bash
            yerr = np.vstack([np.abs(sel - pd.to_numeric(loCI, errors="coerce")),
                              np.abs(pd.to_numeric(upCI, errors="coerce") - sel)])

            fig, ax = plt.subplots(figsize=(11,4))
            sc = ax.scatter(df["POS"], sel, c=-np.log10(df["p_raw"]),
                            vmin=0, vmax=10, s=40, zorder=3)
            ax.errorbar(df["POS"], sel, yerr=yerr,
                        fmt='none', ecolor="grey", elinewidth=1,
                        capsize=1.5, zorder=2)

            texts = []
            for i, r in df.iterrows():
                sig = stars(r["p_raw"])
                if sig:
                    trans_star = offset_copy(ax.transData, fig, y=4, units='points')
                    ax.text(r["POS"], upCI.iloc[i], sig, ha='center', va='bottom',
                            transform=trans_star, color="firebrick", fontsize=7)
                    trans_rs = offset_copy(ax.transData, fig, y=12, units='points')
                    texts.append(ax.text(r["POS"], upCI.iloc[i], r["rsID"],
                                         ha='center', va='bottom',
                                         transform=trans_rs, fontsize=7))

            try:
                adjust_text(texts, arrowprops=dict(arrowstyle="-", color='grey', lw=.5))
            except Exception:
                # se adjust_text dà problemi, continuiamo senza
                pass

            ax.set_xlabel("Genomic position (bp)")
            ax.set_ylabel("Selection coefficient (s)")
            ax.set_title(f"Epoch {k}   (* P<0.05  ** P<0.01  *** P<0.001)")
            ax.grid(alpha=.4)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_: f"{int(x)}"))

            cbar = fig.colorbar(sc, ax=ax, pad=0.01)
            cbar.set_label("−log10 P  (0–10)")

            fn_png = f"{prefix}_epoch{k}.png"
            fn_pdf = f"{prefix}_epoch{k}.pdf"
            fig.tight_layout()
            fig.savefig(fn_png, dpi=300)
            fig.savefig(fn_pdf)
            plt.close(fig)
            print(f"Wrote: {os.path.basename(fn_png)}")

    except ImportError as e:
        warn(f"Plot skipped (missing deps): {e}")

    print("\n")
    info("Phase-2 completed.")
    print()

# ============== CLI ==========================================================

def main():
    p = argparse.ArgumentParser(
        description="CLUES2Companion Phase-2 CLI (Relate/SINGER → CLUES2)"
    )
    p.add_argument("--engine", choices=["Relate","Singer"], required=True,
                   help="Continue with Relate or Singer branch")
    p.add_argument("--chr", required=True, help="Chromosome (e.g. 2, X)")
    p.add_argument("--out_prefix", required=True,
                   help="Population/output prefix used in Phase-1")
    p.add_argument("--phase1_dir", help="Override Phase-1 output directory (auto-detected by default)")

    # CLUES2 params (comuni)
    p.add_argument("--tcutoff", type=int, required=True, help="tCutoff (e.g. 1000)")
    p.add_argument("--df", type=int, required=True, help="df (e.g. 600)")
    p.add_argument("--h", type=float, default=0.5,
               help="Dominance coefficient (default: 0.5, additive model)")
    p.add_argument("--anc_samps", help="AncientSamps file (optional)")
    p.add_argument("--anc_haps",  help="AncientHaps  file (optional)")
    p.add_argument("--noTraj", action="store_true", help="Disable allele trajectory")
    p.add_argument("--CI", type=float, help="Confidence interval (e.g. 0.95)")
    p.add_argument("--timeBins", help="Space/comma separated epoch breakpoints (optional)")


    # Relate-specific
    p.add_argument("--num_samples", type=int, help="Importance sampling of branch lengths (Relate)")
    p.add_argument("--mu", type=float, default=1.25e-8, help="Mutation rate for SampleBranchLengths (Relate)")

    # Singer-specific
    p.add_argument("--region_start", type=int, help="START bp of region (Singer)")
    p.add_argument("--region_end",   type=int, help="END bp of region (Singer)")
    p.add_argument("--Ne", type=int, help="Effective population size for inference.py (Singer)")

    args = p.parse_args()
    run_phase2(args)

if __name__ == "__main__":
    main()


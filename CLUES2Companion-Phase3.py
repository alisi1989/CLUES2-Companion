
#!/usr/bin/env python3
import argparse, os, sys, subprocess, random, glob, json, math, shutil

# ================= Helpers ===================================================

def info(msg): print(f"\033[1;32m→ {msg}\033[0m")
def warn(msg): print(f"\033[33m[WARN] {msg}\033[0m", file=sys.stderr)
def error(msg): sys.exit(f"\033[31m[ERROR] {msg}\033[0m")

def abs_path(p: str) -> str:
    if not p: return ""
    return os.path.abspath(os.path.expanduser(p))

def run_cmd(cmd, cwd=None):
    if isinstance(cmd, (list, tuple)):
        cmd_list = list(map(str, cmd))
    else:
        import shlex
        cmd_list = shlex.split(cmd)
    info(f"Running: {' '.join(cmd_list)}")
    try:
        subprocess.run(cmd_list, check=True, cwd=cwd)
    except subprocess.CalledProcessError:
        error(f"Command failed: {' '.join(cmd_list)}")

def ensure_file(path, msg=None):
    if not os.path.isfile(path):
        error(msg or f"Missing file: {path}")

def ensure_dir(path, msg=None):
    if not os.path.isdir(path):
        error(msg or f"Missing directory: {path}")

def build_breaks(g_start: int, g_end: int, step: int):
    if step <= 0 or g_start >= g_end:
        error("STEP must be > 0 and G_START < G_END")
    breaks = list(range(g_start, g_end + 1, step))
    if breaks[-1] != g_end:
        breaks.append(g_end)
    return breaks

def find_value_cols(header):
    """Trova indici colonne chiave in un'inference.txt (robusto a piccole variazioni)."""
    h2i = {h:i for i,h in enumerate(header)}
    # -log10(p) può chiamarsi "-log10(p-value)" oppure "-log10(p)"
    logp_key = None
    for k in ("-log10(p-value)", "-log10(p)"):
        if k in h2i: logp_key = k; break
    # logLR
    loglr_key = "logLR" if "logLR" in h2i else None
    # EpochX_start/end & SelectionMLEx (preferisci X=2 se presente)
    sel_keys = [ (k, int(k.split("SelectionMLE")[1])) for k in header if k.startswith("SelectionMLE")]
    epoch_start_keys = [ (k, int(k.split("Epoch")[1].split("_")[0])) for k in header if k.startswith("Epoch") and k.endswith("_start")]
    epoch_end_keys   = [ (k, int(k.split("Epoch")[1].split("_")[0])) for k in header if k.startswith("Epoch") and k.endswith("_end")]
    sel_by_k = {knum:kname for (kname, knum) in sel_keys}
    es_by_k  = {knum:kname for (kname, knum) in epoch_start_keys}
    ee_by_k  = {knum:kname for (kname, knum) in epoch_end_keys}
    # Preferisci epoch 2 se completo, altrimenti l'ultimo disponibile
    preferred_k = None
    if 2 in sel_by_k and 2 in es_by_k and 2 in ee_by_k:
        preferred_k = 2
    else:
        avail = sorted(set(sel_by_k) & set(es_by_k) & set(ee_by_k))
        preferred_k = avail[-1] if avail else None
    return dict(
        logp_key=logp_key,
        loglr_key=loglr_key,
        sel_key = sel_by_k.get(preferred_k),
        e_start = es_by_k.get(preferred_k),
        e_end   = ee_by_k.get(preferred_k)
    )

def read_freq_for_rsid(freq_file, rsid):
    """Ritorna la freq (string) per l'rsid (colonna 4 o ultima numerica)."""
    ensure_file(freq_file, f"Frequency file not found: {freq_file}")
    with open(freq_file) as f:
        header = f.readline()
        for ln in f:
            ln = ln.strip()
            if not ln: continue
            parts = ln.split()
            if parts[0] == rsid:
                # usa 4a col se presente, altrimenti l'ultima numerica
                if len(parts) >= 4:
                    return parts[3]
                # fallback: ultima numerica
                for x in reversed(parts):
                    try:
                        float(x); return x
                    except: pass
                break
    error(f"Derived allele frequency for {rsid} not found in {os.path.basename(freq_file)}")

def parse_inference_file(path):
    """Legge un *_inference.txt e ritorna (header:list, data:list)."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        data   = f.readline().rstrip("\n").split("\t")
    return header, data

def last_selection_k(header):
    ks = []
    for h in header:
        if h.startswith("SelectionMLE"):
            try: ks.append(int(h.split("SelectionMLE")[1]))
            except: pass
    return max(ks) if ks else None

def infer_windows_from_prefix(prefix_root):
    """Raccoglie files tipo <prefix>_<L>_<R>_inference.txt → [(L,R,s,pval)]."""
    windows = []
    for path in glob.glob(f"{prefix_root}_*_inference.txt"):
        base = os.path.basename(path)
        # ..._<L>_<R>_inference.txt
        try:
            tail = base.rsplit("_inference.txt", 1)[0]
            left, right = map(int, tail.rsplit("_", 2)[-2:])
        except Exception:
            continue
        header, data = parse_inference_file(path)
        # prendi ultimo SelectionMLE disponibile
        lk = last_selection_k(header)
        if lk is None:
            continue
        sel_name = f"SelectionMLE{lk}"
        try:
            s_val = float(data[header.index(sel_name)])
        except Exception:
            continue
        # pval (se c'è)
        pval = float("nan")
        for key in ("-log10(p-value)", "-log10(p)"):
            if key in header:
                try:
                    logp = float(data[header.index(key)])
                    pval = 10**(-logp)
                except: pass
                break
        windows.append((left, right, s_val, pval))
    # ordina da antico→recente (left decrescente)
    windows.sort(key=lambda t: t[0], reverse=True)
    return windows

def onset_from_windows(windows):
    """Criteri: 4→3→2 finestre consecutive con s>0, poi first positive, poi max s."""
    if not windows:
        return None
    def consecutive(k):
        for i in range(len(windows)-k+1):
            seg = windows[i:i+k]
            if all(w[2] > 0 for w in seg):
                st, en, s, p = seg[0]
                return dict(start=st, end=en, s=s, k=k, method=f"{k} consecutive windows > 0")
        return None
    for k in (4,3,2):
        res = consecutive(k)
        if res: return res
    # fallback: prima positiva
    for st,en,s,p in windows:
        if s > 0:
            return dict(start=st, end=en, s=s, k=1, method="first positive s")
    # fallback: massima s
    st,en,s,p = max(windows, key=lambda t: t[2])
    return dict(start=st, end=en, s=s, k=1, method="max s")

def percentiles(values, q_low=2.5, q_high=97.5):
    try:
        import numpy as np
        arr = np.array(values, dtype=float)
        med = float(np.median(arr))
        lo  = float(np.percentile(arr, q_low))
        hi  = float(np.percentile(arr, q_high))
        return int(round(med)), int(math.floor(lo)), int(math.ceil(hi))
    except Exception:
        # fallback puro-Python
        vals = sorted(map(float, values))
        def q(p):
            idx = (len(vals)-1) * (p/100.0)
            lo = int(math.floor(idx)); hi = int(math.ceil(idx))
            if lo == hi: return vals[lo]
            w = idx - lo
            return vals[lo]*(1-w) + vals[hi]*w
        med = q(50); lo = q(q_low); hi = q(q_high)
        return int(round(med)), int(math.floor(lo)), int(math.ceil(hi))

# ================= Phase-3 Runner (Relate only) ==============================

def run_phase3(args):
    base_dir   = os.getcwd()
    chr_       = args.chr
    out_pref   = args.out_prefix
    rsid       = args.rsid
    df_score   = args.df
    g_start    = args.g_start
    g_end      = args.g_end
    step       = args.step

    print("\n          ╔═════════════════════════════════════════════════════════╗")
    print("          ║        ⏱️  PHASE 3 – Dating a selective sweep  ⏱️         ║")
    print("          ║   Please read the manual carefully before proceeding    ║")
    print("          ╚═════════════════════════════════════════════════════════╝\n")

    # ---- directories (Relate only) ----------------------------------------
    work_base = os.path.join(base_dir, "output_CLUES2Companion-Relate")
    p1_dir = abs_path(args.phase1_dir) if args.phase1_dir else \
             os.path.join(work_base, "phase1", f"{out_pref}_chr{chr_}")
    p2_dir = abs_path(args.phase2_dir) if args.phase2_dir else \
             os.path.join(work_base, "phase2", f"{out_pref}_chr{chr_}")
    ensure_dir(p1_dir, f"Phase-1 output not found: {p1_dir}")
    ensure_dir(p2_dir, f"Phase-2 output not found: {p2_dir}")

    work_dir = os.path.join(work_base, "phase3", f"{out_pref}_chr{chr_}")
    os.makedirs(work_dir, exist_ok=True)
    info(f"Phase 3 working dir: {work_dir}")

    # ---- required Phase-1/2 artifacts -------------------------------------
    coal_file   = os.path.join(p1_dir, f"{out_pref}_EPS4COAL_chr{chr_}.coal")
    ensure_file(coal_file, "*.coal file missing")

    derived_file = os.path.join(p1_dir, f"{out_pref}_Derived_{rsid}.txt")
    ensure_file(derived_file, f"Derived file missing ({derived_file})")

    freq_glob = os.path.join(p1_dir, f"{out_pref}_Frequency_chr{chr_}_*.txt")
    freq_files = sorted(glob.glob(freq_glob))
    if not freq_files: error("frequency file not found in Phase-1 dir")
    freq_file = freq_files[0]
    freq_val  = read_freq_for_rsid(freq_file, rsid)

    times0 = os.path.join(p2_dir, f"{out_pref}_times_chr{chr_}", f"{rsid}_times.txt")
    ensure_file(times0, f"times file missing ({times0})")
    times_file = os.path.join(work_dir, f"{rsid}_times.txt")
    shutil.copyfile(times0, times_file)

    # ---- tool paths --------------------------------------------------------
    relate_root = os.path.join(base_dir, "Relate")
    clues_dir   = os.path.join(base_dir, "CLUES2")
    inf_py      = os.path.join(clues_dir, "inference.py")
    rtc_py      = os.path.join(clues_dir, "RelateToCLUES.py")
    ensure_file(inf_py, "CLUES2/inference.py not found")
    ensure_file(rtc_py, "CLUES2/RelateToCLUES.py not found")

    # ---- output phase-3 (Dating) ------------------------------------------
    dat_dir  = os.path.join(work_dir, "Dating"); os.makedirs(dat_dir, exist_ok=True)
    inf_out_root = os.path.join(dat_dir, rsid)     # prefix for first scan
    json_out1    = os.path.join(dat_dir, f"{rsid}_InitialOnset_Dating.json")

    # ===================== (A) Sliding-window scan ==========================
    breaks = build_breaks(g_start, g_end, step)
    print("\nBreak-points (STEP (left-right)):")
    for i in range(len(breaks)-1):
        print(f"{step:4d} ({breaks[i]}-{breaks[i+1]})")
    print(f"tCutoff      : {breaks[-1]}")
    print(f"Window size  : {step}\n")

    print("▶  Scanning windows:")
    for i in range(len(breaks)-1):
        left  = breaks[i]
        right = breaks[i+1]
        win_pref = f"{inf_out_root}_{left}_{right}"
        info(f"    CLUES2 • window [{left}-{right}]")
        cmd = [
            sys.executable, inf_py,
            "--times",   times_file,
            "--popFreq", str(freq_val),
            "--out",     win_pref,
            "--df",      str(df_score),
            "--CI",      "0.95",
            "--tCutoff", str(right),
            "--coal",    coal_file,
            "--noAlleleTraj",
            "--h",       str(args.h)
        ]
        if left > 0:
            cmd += ["--timeBins", str(left)]
        try:
            run_cmd(cmd)
        except SystemExit:
            warn("      ➜ CLUES2 failed – window skipped")
            continue

    # ---- (A-2) Collect & onset --------------------------------------------
    wins = infer_windows_from_prefix(inf_out_root)
    if not wins:
        error("No inference files produced in the scan")
    onset = onset_from_windows(wins)
    st, en, s_mle = onset["start"], onset["end"], onset["s"]
    med_gen   = int(round((st + en) / 2))
    med_years = med_gen * 28

    print(f"\nInitial onset (median of window) ≈ {med_gen} generations  (≈ {med_years} years)")
    print(f"   window   : {st} – {en} generations")
    print(f"   s(MLE)   : {s_mle:.5f}")
    print(f"   method   : {onset['method']}\n")

    with open(json_out1, "w") as out:
        json.dump(dict(
            rsID               = rsid,
            population         = out_pref,
            chromosome         = chr_,
            median_onset_gen   = med_gen,
            median_onset_years = med_years,
            epoch_start        = st,
            epoch_end          = en,
            epoch_start_years  = st * 28,
            epoch_end_years    = en * 28,
            method             = onset["method"],
            s_MLE              = round(float(s_mle), 5)
        ), out, indent=2)
    info(f"JSON written → {json_out1}")

    # ===================== (B) Bootstrap (optional) =========================
    if not args.bootstrap:
        print("\nBootstrap step skipped – Phase 3 completed.")
        return

    if args.num_samples is None or args.num_samples <= 0:
        error("Bootstrap requires --num_samples (importance samples for SampleBranchLengths)")

    # parametri bootstrap
    b_start = args.boot_start
    b_end   = args.boot_end
    b_step  = args.boot_step
    nboot   = args.nboot
    boot_df = args.boot_df
    mu      = args.mu

    if b_start is None or b_end is None:
        error("Bootstrap requires --boot_start and --boot_end")
    b_breaks = build_breaks(b_start, b_end, b_step)
    print("\nBreak-points (STEP (left-right)):")
    for i in range(len(b_breaks)-1):
        print(f"{b_step:4d} ({b_breaks[i]}-{b_breaks[i+1]})")
    print(f"Replicates   : {nboot}\n")

    boot_dir  = os.path.join(dat_dir, f"bootstrap_{rsid}")
    os.makedirs(boot_dir, exist_ok=True)

    # GS prefix & posizione fisica dell'rs
    gs_prefix = os.path.join(p1_dir, f"{out_pref}_GS+COAL_chr{chr_}")
    snp_glob  = os.path.join(p1_dir, f"{out_pref}_SNPs_chr{chr_}_*.txt")
    snp_files = sorted(glob.glob(snp_glob))
    if not snp_files: error("SNPs file not found in Phase-1 dir")
    # trova POS dell'rs
    rs_pos = None
    for sf in snp_files:
        with open(sf) as f:
            for ln in f:
                ln = ln.strip()
                if not ln: continue
                parts = ln.split()
                if parts[0] == rsid:
                    rs_pos = parts[1]
                    break
        if rs_pos: break
    if rs_pos is None:
        error(f"bp position for {rsid} not found in SNPs files")

    # tool script per SampleBranchLengths e RelateToCLUES
    sbl_sh = os.path.join(relate_root, "scripts", "SampleBranchLengths", "SampleBranchLengths.sh")
    if not os.path.exists(sbl_sh):
        error(f"Relate SampleBranchLengths script not found: {sbl_sh}")

    print(f"\n▶  Generating {nboot} bootstrap trees & CLUES2 runs")
    for rep in range(1, nboot+1):
        seed = random.randint(1, 2_000_000_000)
        info(f"[bootstrap {rep}]  SampleBranchLengths (seed={seed})")
        nw_pref = os.path.join(boot_dir, f"nw_{rep}")
        # SampleBranchLengths
        cmd_sbl = [
            sbl_sh,
            "--input",  gs_prefix,
            "--output", nw_pref,
            "--first_bp", rs_pos,
            "--last_bp",  rs_pos,
            "--format", "n",
            "--num_samples", str(args.num_samples),
            "--coal",  coal_file,
            "--seed",  str(seed),
            "-m",      str(mu),
            
        ]
        try:
            run_cmd(cmd_sbl, cwd=boot_dir)
        except SystemExit:
            warn(f"[bootstrap {rep}] SampleBranchLengths failed – skipping replicate")
            continue

        # RelateToCLUES (times)
        info("  RelateToCLUES")
        cmd_rtc = [sys.executable, rtc_py,
                   "--DerivedFile",  derived_file,
                   "--RelateSamples", f"{nw_pref}.newick",
                   "--out", nw_pref]
        try:
            run_cmd(cmd_rtc)
        except SystemExit:
            warn(f"[bootstrap {rep}] RelateToCLUES failed – skipping replicate")
            continue

        # CLUES2 sulle finestre bootstrap
        for i in range(len(b_breaks)-1):
            left, right = b_breaks[i], b_breaks[i+1]
            win_pref = f"{nw_pref}_{left}_{right}"
            info(f"    CLUES2 inference • window [{left} - {right}]")
            cmd_inf = [
                sys.executable, inf_py,
                "--times", f"{nw_pref}_times.txt",
                "--popFreq", str(freq_val),
                "--out", win_pref,
                "--df",  str(boot_df),
                "--CI",  "0.95",
                "--tCutoff", str(right),
                "--coal",    coal_file,
                "--noAlleleTraj",
                "--h",       str(args.h)
            ]
            if left > 0:
                cmd_inf += ["--timeBins", str(left)]
            try:
                run_cmd(cmd_inf)
            except SystemExit:
                warn("      ➜ CLUES2 failed – window skipped")
                continue

        # Dopo tutte le finestre: concat Epoch-2 (o ultimo epoch disponibile)
        combo = os.path.join(boot_dir, f"bootstrap_{rep}_{rsid}.txt")
        with open(combo, "w") as out:
            out.write("logLR\t-log10(p-value)\tEpoch2_start\tEpoch2_end\tSelectionMLE2\n")
            for path in glob.glob(f"{nw_pref}_*_inference.txt"):
                try:
                    header, data = parse_inference_file(path)
                except Exception:
                    continue
                cols = find_value_cols(header)
                # Se manca Epoch2/SelectionMLE2, cade sul'ultimo epoch disponibile
                sel_k  = cols["sel_key"]
                e_s    = cols["e_start"]
                e_e    = cols["e_end"]
                logp_k = cols["logp_key"] or "-log10(p-value)"  # intestazione di output
                loglr_k= cols["loglr_key"]
                if not (sel_k and e_s and e_e and logp_k and loglr_k):
                    continue
                try:
                    row = [
                        data[header.index(loglr_k)],
                        data[header.index(cols["logp_key"])] if cols["logp_key"] in header else "NA",
                        data[header.index(e_s)],
                        data[header.index(e_e)],
                        data[header.index(sel_k)]
                    ]
                except Exception:
                    continue
                out.write("\t".join(row) + "\n")
                # pulizia file finestra per risparmiare spazio
                try: os.remove(path)
                except: pass

    # ===================== (D) Onset per ogni bootstrap ======================
    print("\n▶  Computing onset for each bootstrap…")
    onset_vals = []
    method_vals = []
    for rep in range(1, nboot+1):
        combo = os.path.join(boot_dir, f"bootstrap_{rep}_{rsid}.txt")
        if not os.path.isfile(combo) or os.path.getsize(combo) == 0:
            warn(f"Bootstrap {rep} – file {os.path.basename(combo)} empty or missing, skip")
            continue
        # Leggi finestre, ordina antico→recente
        wins = []
        with open(combo) as f:
            header = f.readline().rstrip("\n").split("\t")
            h2i = {h:i for i,h in enumerate(header)}
            need = ["Epoch2_start","Epoch2_end","SelectionMLE2"]
            if not all(n in h2i for n in need):
                # se header non ha 2, prova a dedurre epoch preferito dal combo (ultimo colonna che inizia per SelectionMLE)
                # ma dato che combo è fissato, se non combacia salta
                warn(f"Bootstrap {rep} – unexpected columns in combo, skip")
                continue
            for ln in f:
                ln = ln.strip()
                if not ln: continue
                parts = ln.split("\t")
                try:
                    st = int(float(parts[h2i["Epoch2_start"]]))
                    en = int(float(parts[h2i["Epoch2_end"]]))
                    s  = float(parts[h2i["SelectionMLE2"]])
                    wins.append((st,en,s))
                except:
                    continue
        wins.sort(key=lambda t: t[0], reverse=True)

        # onset: 4→3→2, poi first positive, poi max
        def med(st,en): return int(round((st+en)/2))
        def find_k(k):
            for i in range(len(wins)-k+1):
                seg = wins[i:i+k]
                if all(w[2] > 0 for w in seg):
                    st,en,_ = seg[0]
                    return med(st,en), f"{k} consecutive s>0"
            return None, None
        onset, method = None, None
        for k in (4,3,2):
            onset, method = find_k(k)
            if onset is not None: break
        if onset is None:
            for st,en,s in wins:
                if s>0:
                    onset, method = med(st,en), "first positive s"
                    break
        if onset is None and wins:
            st,en,_ = max(wins, key=lambda t: t[2])
            onset, method = med(st,en), "max s"
        if onset is None:
            warn(f"Bootstrap {rep} – no onset inferred, skip")
            continue

        onset_vals.append(onset); method_vals.append(method)
        line = f"Bootstrap {rep:2d} : {onset} gen  ({method})"
        print(line)
        # append to a log file in work_dir
        with open(os.path.join(work_dir, "onset_bootstraps.txt"), "a") as lg:
            lg.write(line+"\n")

    # ===================== (E) Statistiche & JSON ============================
    print("\n▶  Summary statistics")
    print(f"Valid onsets collected: {len(onset_vals)}")
    if len(onset_vals) < 2:
        error("Fewer than 2 valid onsets – CI cannot be computed")

    MEDIAN, CI_LOW, CI_HIGH = percentiles(onset_vals, 2.5, 97.5)
    print(f"95 % CI      : {CI_LOW} – {CI_HIGH} gen "
          f"(~{CI_LOW*28} – {CI_HIGH*28} yr)")

    # metodo prevalente
    prev_method = "NA"
    if method_vals:
        from collections import Counter
        c = Counter(method_vals)
        m, n = c.most_common(1)[0]
        prev_method = f"{m} ({n}×)"

    out_json = os.path.join(dat_dir, f"{rsid}_Boostraps_onset_Dating+CI.json")
    with open(out_json, "w") as out:
        json.dump(dict(
            rsID           = rsid,
            population     = out_pref,
            chromosome     = chr_,
            CI95_low_gen   = int(CI_LOW),
            CI95_high_gen  = int(CI_HIGH),
            CI95_low_year  = int(CI_LOW)*28,
            CI95_high_year = int(CI_HIGH)*28,
            bootstraps     = int(nboot),
            step           = int(b_step),
            method         = prev_method
        ), out, indent=2)
    print(f"JSON saved → {out_json}")
    print("\nPhase-3 completed.")

# ================= CLI =======================================================

def main():
    p = argparse.ArgumentParser(
        description="CLUES2Companion Phase-3 CLI (Relate only) – Sliding-window dating + optional bootstrap"
    )
    p.add_argument("--chr", required=True, help="Chromosome (e.g. 2, X)")
    p.add_argument("--out_prefix", required=True, help="Population/output prefix used in Phase-1/2")
    p.add_argument("--rsid", required=True, help="Target SNP rsID to date (e.g. rs123)")

    # Scan (A)
    p.add_argument("--df", type=int, required=True, help="df score for CLUES2 (e.g. 600)")
    p.add_argument("--g_start", type=int, required=True, help="Initial epoch to scan (generations ago)")
    p.add_argument("--g_end",   type=int, required=True, help="Final epoch to scan (generations ago)")
    p.add_argument("--step",    type=int, default=50,    help="Window size (non-overlapping), default 50")

    # Paths override (optional)
    p.add_argument("--phase1_dir", help="Override Phase-1 directory (auto-detected by default)")
    p.add_argument("--phase2_dir", help="Override Phase-2 directory (auto-detected by default)")

    # Bootstrap (B)
    p.add_argument("--bootstrap", action="store_true", help="Enable bootstrap dating")
    p.add_argument("--boot_start", type=int, help="Bootstrap scan START (generations ago)")
    p.add_argument("--boot_end",   type=int, help="Bootstrap scan END (generations ago)")
    p.add_argument("--boot_step",  type=int, default=25, help="Bootstrap window size (default 25)")
    p.add_argument("--nboot",      type=int, default=100, help="Number of bootstrap replicates (default 100)")
    p.add_argument("--boot_df",    type=int, default=450, help="df score for CLUES2 in bootstrap (default 450)")
    p.add_argument("--num_samples", type=int, help="Importance sampling for SampleBranchLengths (required if --bootstrap)")
    p.add_argument("--mu", type=float, default=1.25e-8, help="Mutation rate for SampleBranchLengths (default 1.25e-8)")
    p.add_argument("--h", type=float, default=0.5,
               help="Dominance coefficient (default: 0.5, additive model)")


    args = p.parse_args()
    run_phase3(args)

if __name__ == "__main__":
    main()

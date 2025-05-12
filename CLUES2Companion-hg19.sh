#!/usr/bin/env bash

# ------------------------------------------------------------------
#  CLUES2Companion pipeline – v2025‑04‑16
#  • Phase‑1  : Relate + SNP extraction + Derived files + frequency
#  • Phase‑2  : BranchLengths → RelateToCLUES → inference → merge
#  • Phase‑3  : Dating target SNP(s)
# ------------------------------------------------------------------
set -euo pipefail
IFS=$'\n\t'

# ---------- helpers ----------------------------------------------------------
abs_path() {
  local p="$1"
  [[ -z "$p" ]] && { echo ""; return; }
  case "$p" in
    /*) echo "$p" ;;
    ~*) echo "${p/#\~/${HOME}}" ;;
     *) echo "$(pwd)/$p" ;;
  esac
}

SCRIPT_NAME=$(basename "$0")
BASE_DIR=$(pwd)
WORK_BASE="${BASE_DIR}/output_C2Companion"; mkdir -p "$WORK_BASE"
LOG_DIR="${WORK_BASE}/log";          mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/${SCRIPT_NAME%.*}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# colours
if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
  BOLD="\e[1m"; GREEN="\e[32m"; YELLOW="\e[33m"; RED="\e[31m"; RESET="\e[0m"
else
  BOLD=""; GREEN=""; YELLOW=""; RED=""; RESET=""
fi
info()  { echo -e "${BOLD}${GREEN}→ $*${RESET}"; }
warn()  { echo -e "${YELLOW}[WARN] $*${RESET}"; }
error() { echo -e "${RED}[ERROR] $*${RESET}"; exit 1; }

mark_done(){ touch "$1/.done"; }
phase_done(){ [[ -f "$1/.done" ]]; }

# ---------- menu -------------------------------------------------------------
echo
echo "******  CLUES2Companion – please cite CLUES2 and CLUES2Companion  ******"
cat <<EOF
Choose phase to run
  1) Phase-1  : Relate (*.mut, *.anc, *.coal files) and SNP/Derived/DAF
  2) Phase-2  : Relate (BranchLengths) → RelateToCLUES.py → inference.py → outputs
  3) Phase-3  : Dating target SNP(s)
EOF
read -rp "Enter option (1/2/3): " OPTION
echo

# ---------- PHASE‑1  ---------------------------------------------------------
phase1() {
  PH="phase1"
  echo -e "\n          ╔═════════════════════════════════════════════════════════╗"
  echo -e "          ║         🚀  PHASE 1: RELATE & SNP EXTRACTION  🚀        ║"
  echo -e "          ║   Please read the manual carefully before proceeding    ║"
  echo -e "          ╚═════════════════════════════════════════════════════════╝\n"

  # 1. User input
  # -------------------------------------------------------------------------
  read -rp  "Choose the chromosome to analyze (e.g. 2, 17, X): " chr
  read -e -rp "Prefix of phased vcf/bcf file (e.g. example/Finnish without _chrN.vcf.gz): " vcf_prefix_in
  read -e -rp "Path to population‑labels file (*.poplabels)(e.g. example/Finnish.poplabels): "  pop_labels_in
  read -rp  "Start bp of target region: " start_bp
  read -rp  "End bp of target region: " end_bp
  read -e -rp "Prefix of output name (e.g. Finnish): " out_prefix
  export OUT_PREFIX="$out_prefix"

  RUN_ID="${out_prefix}_chr${chr}"
  WORK_DIR="${WORK_BASE}/${PH}/${RUN_ID}"
  mkdir -p "$WORK_DIR"
  echo
  echo -e "[INFO] Working dir: $WORK_DIR"

 
  # Relate path is fixed
  relate_root="${BASE_DIR}/Relate"
  [[ -x "${relate_root}/bin/Relate" ]] || error "Relate binary not found in ${relate_root}/bin"

  # -------------------------------------------------------------------------
  # 2. required_files files
  # -------------------------------------------------------------------------
  required_files_dir="${BASE_DIR}/required_files-hg19"
  ancestral="${required_files_dir}/ancestor/human_ancestor_chr${chr}.fa"
  mask="${required_files_dir}/mask/PilotMask_chr${chr}.fasta"
  genmap="${required_files_dir}/map/genetic_map_chr${chr}.txt"
  for f in "$ancestral" "$mask" "$genmap"; do [[ -f $f ]] || error "Missing $f"; done
  pop_labels=$(abs_path "$pop_labels_in")

  # -------------------------------------------------------------------------
  # 3. internal prefixes
  # -------------------------------------------------------------------------
  vcf_gz="${vcf_prefix_in}_chr${chr}"             # senza .vcf.gz
  haps_path="${WORK_DIR}/${out_prefix}_chr${chr}.haps"
  sample_path="${WORK_DIR}/${out_prefix}_chr${chr}.sample"
  pif_base="${WORK_DIR}/${out_prefix}_PIF_chr${chr}"
  gs_base="${WORK_DIR}/${out_prefix}_GS_chr${chr}"
  gs_coal="${WORK_DIR}/${out_prefix}_GS+COAL_chr${chr}"
  coal_base="${WORK_DIR}/${out_prefix}_EPS4COAL_chr${chr}"
  
  # spinner initializing#
  spinner_bar() {
  local pid=$1
  local delay=${2:-1}
  local max=${3:-30}
  local count=0
  echo "THIS MAY TAKE A WHILE, PLEASE WAIT."
  printf "["

  # finché il processo esiste…
  while kill -0 "$pid" 2>/dev/null; do
    (( count++ ))
    if (( count > max )); then
      printf "\r[\033[K"
      count=1
    fi
    printf "▊"
    sleep "$delay"
  done
  # pulisci tutta la riga e vai a capo
  printf "] done!"
}


  
  # -------------------------------------------------------------------------
  # 4. Convert VCF
  # -------------------------------------------------------------------------
  echo
  info "Convert vcf → *.haps and *.sample files"
  "${relate_root}/bin/RelateFileFormats" --mode ConvertFromVcf \
      -i "$vcf_gz" --haps "$haps_path" --sample "$sample_path" -o "${WORK_DIR}/${out_prefix}_chr${chr}" 2>/dev/null
  echo "✓ Convert vcf file done"
  # -------------------------------------------------------------------------
  # 5. PrepareInputFiles
  # -------------------------------------------------------------------------
  
  echo
  info "Preparing input files → *.haps.gz, *.sample.gz, *.dist and *.annot"
  
  "${relate_root}/scripts/PrepareInputFiles/PrepareInputFiles.sh" \
      --haps      "$haps_path" \
      --sample    "$sample_path" \
      --ancestor  "$ancestral" \
      --mask      "$mask" \
      --poplabels "$pop_labels" \
      -o          "$pif_base" 2>/dev/null &
  
  pid=$!
  spinner_bar "$pid" 0.5 30
  wait "$pid" || true
  echo "✓ Preparing input files done"
  
  pif_haps="${pif_base}.haps.gz"
  pif_sample="${pif_base}.sample.gz"
  pif_annot="${pif_base}.annot"

# -------------------------------------------------------------------------
# 6. Relate run-1  (random Ne)
# -------------------------------------------------------------------------
echo
info "Running Relate mode All with random Ne → *.anc and *.mut files"

pushd "$WORK_DIR" >/dev/null   # entra nella dir di lavoro

"${relate_root}/bin/Relate" --mode All \
    --haps  "$pif_haps" --sample "$pif_sample" \
    --map   "$genmap"   --annot "$pif_annot" \
    -N 30000 -m 1.25e-8 --seed 1 -o GS_tmp 2>/dev/null &

  pid=$!
  spinner_bar "$pid" 0.5 30
  wait "$pid" || true

  for f in "$WORK_DIR"/GS_tmp.*; do mv "$f" "${gs_base}_run1${f##*/GS_tmp}"; done
  echo -e "✓ Relate mode All with random Ne done"
  popd >/dev/null


  # -------------------------------------------------------------------------
  # 7. EstimatePopulationSize  → .coal
  # -------------------------------------------------------------------------
  echo
  info "Estimating Population size using Relate → *.coal file"
  
  pushd "$WORK_DIR" >/dev/null
  
  "${relate_root}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh" \
        -i "${gs_base}_run1" -o "${coal_base}" \
        --noanc 1 --poplabels "$pop_labels" \
        -m 1.25e-8 --years_per_gen 28 --seed 1 2>/dev/null &
      

  pid=$!
  spinner_bar "$pid" 0.5 30
  wait "$pid" || true
 if ! wait    "$pid"; then
   warn "plot failed (missing R/ggplot); skipping plot and continuing."
 fi
 
 popd >/dev/null
 echo "✓ Estimating Population size done"
 
  # -------------------------------------------------------------------------
  # 8. Relate run‑2  (final .anc/.mut)
  # -------------------------------------------------------------------------
 echo
 info "Running Relate mode All with *.coal file → *.anc and *.mut based on coalescence"
  
 pushd "$WORK_DIR" >/dev/null   # entra nella dir di lavoro

"${relate_root}/bin/Relate" --mode All \
    --haps  "$pif_haps" --sample "$pif_sample" \
    --map   "$genmap"   --annot "$pif_annot" \
    --coal "${coal_base}.coal" \
    -m 1.25e-8 --seed 1 -o GS_tmp2 2>/dev/null &
 
  pid=$!
  spinner_bar "$pid" 0.5 30
  wait "$pid" || true
  for f in GS_tmp2.*; do mv "$f" "${gs_coal}${f#GS_tmp2}"; done
  popd >/dev/null
  echo "✓ Relate mode All with *.coal file done"
  
  mut_file="${gs_coal}.mut"
  haps_file="${pif_base}.haps.gz"

  # -------------------------------------------------------------------------
  # 9. SNP extraction
  # -------------------------------------------------------------------------
  echo
  info "Extracting SNPs using cyvcf2"
  export OUTPUT_DIR="$WORK_DIR" CHR="$chr" VCF_GZ="${vcf_gz}.vcf.gz"
  export START_BP="$start_bp" END_BP="$end_bp"
  snp_file="${WORK_DIR}/${OUT_PREFIX}_SNPs_chr${chr}_${start_bp}_${end_bp}.txt"
  export SNP_FILE="$snp_file"

  python - <<'PY'
import os
from cyvcf2 import VCF
vcf   = VCF(os.environ["VCF_GZ"])
reg   = f"chr{os.environ['CHR']}:{os.environ['START_BP']}-{os.environ['END_BP']}"
outp  = os.environ["SNP_FILE"]; cnt=0
with open(outp,"w") as out:
    for v in vcf(reg):
        if v.ID and v.ID!=".":
            out.write(f"{v.ID}\t{v.POS}\n"); cnt+=1
print(f"Extracted {cnt} SNPs → {outp}")
PY
  [[ -s "$snp_file" ]] || error "No SNPs extracted."

  # -------------------------------------------------------------------------
  # 10. Derived alleles & frequency
  # -------------------------------------------------------------------------
  echo
  info "Polarizing derived alleles and calculating derived allele frequency"
  export HAPS_FILE="$haps_file" MUT_FILE="$mut_file"
  freq_file="${WORK_DIR}/${OUT_PREFIX}_Frequency_chr${chr}_${start_bp}_${end_bp}.txt"
  export FREQ_FILE="$freq_file"

  python - <<'PY'
import os, gzip
op = lambda p,m="rt": gzip.open(p,m) if p.endswith(".gz") else open(p,m)
prefix=os.environ["OUT_PREFIX"]; outdir=os.environ["OUTPUT_DIR"]

# snp list
snps=[ln.split() for ln in open(os.environ["SNP_FILE"])]
# mut dict
mut={}
with op(os.environ["MUT_FILE"]) as f:
    next(f)
    for ln in f:
        fs=ln.strip().split(";")
        if len(fs)<11 or fs[3]==".":continue
        rs=fs[3]; a,d=fs[10].split("/")
        if fs[7]=="1": a,d=d,a
        mut[rs]=d.upper()
# haps dict
haps={}
with op(os.environ["HAPS_FILE"]) as f:
    for ln in f:
        p=ln.split();
        if len(p)<6: continue
        haps[p[1]]=(p[3].upper(),p[4].upper(),p[5:])

out=open(os.environ["FREQ_FILE"],"w")
out.write("rsID\tPOS\tallele_der\tfreq_der\n")
for rs,pos in snps:
    if rs not in mut or rs not in haps: continue
    der=mut[rs]; ref,alt,vec=haps[rs]
    bits=[int(b) if alt==der else 1-int(b) for b in vec]
    freq=sum(bits)/len(bits)
    dpath=f"{outdir}/{prefix}_Derived_{rs}.txt"
    with open(dpath,"w") as d: d.write("\n".join(map(str,bits)))
    out.write(f"{rs}\t{pos}\t{der}\t{freq:.4f}\n")
    print(f"  ↳ {dpath}")
print(f"✓ Frequency table → {os.environ['FREQ_FILE']}")
PY

# --- clean‑up: remove intermediate files we no longer need -------------------
rm -f "${WORK_DIR}"/*.rate \
      "${WORK_DIR}"/*.pairwise.coal \
      "${WORK_DIR}"/*.pairwise.bin \
      "${WORK_DIR}"/*.annot \
      "${WORK_DIR}"/*.haps \
      "${WORK_DIR}"/*.sample \
      "${WORK_DIR}"/*_run1.anc \
      "${WORK_DIR}"/*_run1.mut
      
  mark_done "$WORK_DIR"
  echo
  info "✓ Phase‑1 completed — outputs in ${WORK_DIR}"
  printf "\n"
}

# ---------- PHASE-2  (BranchLengths → RelateToCLUES → CLUES → merge) ----------
phase2() {

  PH="phase2"

  echo -e "\n          ╔═════════════════════════════════════════════════════════╗"
  echo -e "          ║        🧬  PHASE-2 – (requires Phase‑1 outputs)  🧬     ║"
  echo -e "          ║   Please read the manual carefully before proceeding    ║"
  echo -e "          ╚═════════════════════════════════════════════════════════╝\n"
  # -------------------------------------------------------------------------
  # 1.  INPUTS
  # -------------------------------------------------------------------------
  read -rp  "Choose the chromosome to analyze (e.g. 2, 17, X): " chr
  read -e -rp "Enter population prefix used in Phase‑1 (e.g. Finnish): " out_prefix
  export OUT_PREFIX="$out_prefix"
  
  RUN_ID="${OUT_PREFIX}_chr${chr}"
  WORK_DIR="${WORK_BASE}/${PH}/${RUN_ID}"
  mkdir -p "$WORK_DIR"

  # where to find Phase-1 output (default location is output_C2Companion/phase1)
  phase1_dir_default="${WORK_BASE}/phase1/${out_prefix}_chr${chr}"
  read -e -rp "Phase‑1 auto-detect directory [ENTER = ${phase1_dir_default}] or provide a different folder with phase1 outputs :" phase1_dir
  phase1_dir=${phase1_dir:-$phase1_dir_default}
  [[ -d "$phase1_dir" ]] || error "Phase‑1 directory not found: $phase1_dir"

  # -------------------------------------------------------------------------
  # 2. finding frequency and SNPs files
  # -------------------------------------------------------------------------
  freq_file=$(ls "${phase1_dir}/${OUT_PREFIX}_Frequency_chr${chr}_"*.txt 2>/dev/null | head -1)
  [[ -f "$freq_file" ]] || error "Frequency file not found in ${phase1_dir}"
  echo
  info "Using frequency file: $freq_file"

  snp_file=$(ls "${phase1_dir}/${OUT_PREFIX}_SNPs_chr${chr}_"*.txt 2>/dev/null | head -1)
  [[ -f "$snp_file" ]] || error "SNP coordinate file not found in ${phase1_dir}"
  info "Using SNPs file: $snp_file"
  echo
  # -------------------------------------------------------------------------
  # 3.  Relate (.anc/.mut) and .coal done by Phase1
  # -------------------------------------------------------------------------
  gs_prefix="${phase1_dir}/${OUT_PREFIX}_GS+COAL_chr${chr}"
  [[ -f "${gs_prefix}.anc" || -f "${gs_prefix}.anc.gz" ]] \
        || error "*.anc file not found for prefix $gs_prefix"

  coal_file="${phase1_dir}/${OUT_PREFIX}_EPS4COAL_chr${chr}.coal"
  [[ -f "$coal_file" ]] || error "*.coal file not found: $coal_file"

  # -------------------------------------------------------------------------
  # 4.  Derived files glob
  # -------------------------------------------------------------------------
  derived_glob="${phase1_dir}/${OUT_PREFIX}_Derived_*.txt"
  compgen -G "$derived_glob" >/dev/null || error "Derived files not found in ${phase1_dir}"

  # -------------------------------------------------------------------------
  # 5.  third part scripts sources
  # -------------------------------------------------------------------------
  relate_root="$(dirname "$0")/Relate"
  [[ -x "${relate_root}/bin/Relate" ]] || error "Relate binaries not found"

  # CLUES‑v2 helper scripts sono sempre in  ./CLUES2/
  clues_dir="$(dirname "$0")/CLUES2"
  rtc_py="${clues_dir}/RelateToCLUES.py"
  inf_py="${clues_dir}/inference.py"

  [[ -f "$rtc_py" ]] || error "RelateToCLUES.py not found in $clues_dir"
  [[ -f "$inf_py" ]] || error "inference.py not found in $clues_dir"

  # -------------------------------------------------------------------------
  # 6.  CLUES2 parameters
  # -------------------------------------------------------------------------
  read -rp "tCutoff (e.g. 1000): " tcutoff
  read -rp "df (e.g. 600): "      df_score
  read -rp "importance sampling of branch lengths: " num_samples
  read -rp "AncientSamps file (optional, ENTER to skip): " anc_samps
  read -rp "AncientHaps  file (optional, ENTER to skip): " anc_haps
  read -rp "Disable allele trajectory? (y/N): " no_traj
  read -rp "Confidence interval (e.g. 0.95) [ENTER to skip]: " ci_val
  read -rp "TimeBins (a list of epoch breakpoints; optional, ENTER to skip): " time_bins

 #############################################################################
# 7. Directory output Phase‑2
#############################################################################
TREE_DIR="${WORK_DIR}/${OUT_PREFIX}_trees_chr${chr}"
TIMES_DIR="${WORK_DIR}/${OUT_PREFIX}_times_chr${chr}"
INFER_DIR="${WORK_DIR}/${OUT_PREFIX}_inference_chr${chr}"
mkdir -p "$TREE_DIR" "$TIMES_DIR" "$INFER_DIR"

# ------------ spinner (resta invariato) ------------------------------------
spinner_bar() {
  local pid=$1 delay=${2:-1} max=${3:-30} count=0
  echo "THIS MAY TAKE A WHILE, PLEASE WAIT."
  printf "["
  while kill -0 "$pid" 2>/dev/null; do
    (( ++count > max )) && { printf "\r[\033[K"; count=1; }
    printf "▊"; sleep "$delay"
  done
  printf "] done!\n"
}

# ------------ helper: esegui un blocco + spinner + warning -----------------
run_step() {                       # $1 = nome‑funzione, $2/$3 = delay/max
  local func=$1 delay=${2:-0.5} max=${3:-30}
  local tmpwarn; tmpwarn=$(mktemp)

  (
    #  ridefinisco warn solo dentro al subshell
    warn() { echo "[WARN] $*" >>"$tmpwarn"; }
    "$func"
  ) &
  local step_pid=$!
  spinner_bar "$step_pid" "$delay" "$max"
  wait "$step_pid" || true
  [[ -s $tmpwarn ]] && cat "$tmpwarn"
  rm -f "$tmpwarn"
}

#############################################################################
# --- Step A – SampleBranchLengths ------------------------------------------
#############################################################################
stepA() {
  grep -v -E '^(rsID|rsid|#)' "$snp_file" |
  while read -r rsid pos; do
    [[ -z $rsid || -z $pos ]] && continue
    "${relate_root}/scripts/SampleBranchLengths/SampleBranchLengths.sh" \
        --input "$gs_prefix" \
        --output "${TREE_DIR}/${rsid}" \
        --first_bp "$pos" --last_bp "$pos" \
        --format n --num_samples "$num_samples" \
        --coal "$coal_file" --seed 1 -m 1.25e-8 &>/dev/null \
        || warn "BranchLengths failed for $rsid"
  done
}

echo; info "Step A: SampleBranchLengths"
run_step stepA 0.5 30
echo "✓ Step A completed"

#############################################################################
# --- Step B – RelateToCLUES -----------------------------------------------
#############################################################################
stepB() {
  for dfile in $derived_glob; do
    rsid=$(basename "$dfile" | sed "s/^${OUT_PREFIX}_Derived_//; s/\.txt$//")
    nwx="${TREE_DIR}/${rsid}.newick"
    if [[ -f $nwx ]]; then
      python "$rtc_py" --DerivedFile "$dfile" --RelateSamples "$nwx" \
                       --out "${TIMES_DIR}/${rsid}" &>/dev/null \
                       || warn "RTC failed for $rsid"
    else
      warn "newick not found for $rsid"
    fi
  done
}

echo; info "Step B: RelateToCLUES"
run_step stepB 0.5 30
echo "✓ Step B completed"

#############################################################################
# --- Step C – inference.py -------------------------------------------------
#############################################################################
stepC() {
  grep -v -E '^(rsID|rsid|#)' "$freq_file" |
  while read -r rs pos _ freq; do
    tfile="${TIMES_DIR}/${rs}_times.txt"
    if [[ -f $tfile ]]; then
      cmd=( python "$inf_py" --times "$tfile" --popFreq "$freq" \
            --out "${INFER_DIR}/${rs}" --tCutoff "$tcutoff" --df "$df_score" )
      [[ -f $coal_file        ]] && cmd+=( --coal "$coal_file" )
      [[ -n ${ci_val:-}       ]] && cmd+=( --CI "$ci_val" )
      [[ -f ${anc_samps:-}    ]] && cmd+=( --ancientSamps "$anc_samps" )
      [[ -n ${time_bins:-}    ]] && cmd+=( --timeBins "$time_bins" )
      [[ -f ${anc_haps:-}     ]] && cmd+=( --ancientHaps "$anc_haps" )
      [[ $no_traj =~ ^[yY]$   ]] && cmd+=( --noAlleleTraj )
      "${cmd[@]}" &>/dev/null || warn "inference failed for $rs"
    else
      warn "times not found for $rs"
    fi
  done
}

echo; info "Step C: inference.py"
run_step stepC 0.5 30
echo "✓ Step C completed"


  # -------------------------------------------------------------------------
  # Step D – merge finale
  # -------------------------------------------------------------------------
  echo
  info "Step D: merge inference outputs"
  merged="${WORK_DIR}/${OUT_PREFIX}_merged_inference_chr${chr}.tsv"

  # header
  header="rsID\tPOS\tder_freq\tlogLR\t-log10(p)\tEpoch1_start\tEpoch1_end\tSelectionMLE1"
  if [[ -n "${ci_val:-}" ]]; then
      first_ci=$(ls "${INFER_DIR}"/*_CI.txt 2>/dev/null | head -1)
      [[ -n "$first_ci" ]] && header+="\t$(head -1 "$first_ci" | cut -f2-)"
  fi
  echo -e "$header" > "$merged"

  grep -v -E '^(rsID|rsid|#)' "$freq_file" |
  while read -r rs pos _ freq; do
      inf_file="${INFER_DIR}/${rs}_inference.txt"
      [[ -f "$inf_file" ]] || { warn "inference missing for $rs"; continue; }
      core=$(sed -n '2p' "$inf_file")
      ci_part=""
      if [[ -n "${ci_val:-}" ]]; then
          ci_file="${INFER_DIR}/${rs}_CI.txt"
          [[ -f "$ci_file" ]] && ci_part=$(sed -n '2p' "$ci_file" | cut -f2-)
      fi
      echo -e "${rs}\t${pos}\t${freq}\t${core}\t${ci_part}" >> "$merged"
  done
  
  echo
  info "Merged file → $merged"
  mark_done "$WORK_DIR"
  echo
  info "Generating integrated plot"
  export MERGED_TSV="$merged"
  export PLOT_PREFIX="${WORK_DIR}/${OUT_PREFIX}_clues"
  
    python - <<'PYCODE'
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.transforms import offset_copy
from adjustText import adjust_text
from matplotlib.ticker import FuncFormatter

# --- Bash parameters ---
merged_tsv = os.environ["MERGED_TSV"]
prefix     = os.environ["PLOT_PREFIX"]

# --- Useful functions ---
def bh_fdr(pvals):
    p = np.asarray(pvals)
    n = len(p)
    order = p.argsort()
    ranked = np.empty(n, dtype=float)
    ranked[order] = np.minimum.accumulate((p[order] * n) / (np.arange(n,0,-1)))
    return np.clip(ranked, 0, 1)

def stars(p):
    return "***" if p<0.001 else ("**" if p<0.01 else ("*" if p<0.05 else ""))

# --- Reading data ---
df = pd.read_csv(merged_tsv, sep="\t")
df["p_raw"] = 10**(-df["-log10(p)"])
df["p_use"] = df["p_raw"]
df["sig"]   = df["p_use"].apply(stars)

# Finding the "_lower" bound of confidence interval
lower_cols = [c for c in df.columns if c.endswith("_lower")]
if len(lower_cols) != 1:
    raise ValueError(f"Found {len(lower_cols)} '_lower' cols, can't pick one.")
lower_col = lower_cols[0]
upper_col = lower_col.replace("_lower", "_upper")

# Dynamic calculation for yerr
lower_err = np.abs(df["SelectionMLE1"] - df[lower_col])
upper_err = np.abs(df[upper_col]      - df["SelectionMLE1"])
yerr = np.vstack([lower_err.fillna(0), upper_err.fillna(0)])

# --- Drawing ---
fig, ax = plt.subplots(figsize=(11,5))
sc = ax.scatter(df["POS"], df["SelectionMLE1"],
                c=-np.log10(df["p_use"]), cmap="Reds", s=40,
                vmin=0, vmax=10, zorder=3)
ax.errorbar(df["POS"], df["SelectionMLE1"], yerr=yerr,
            fmt='none', ecolor="grey", elinewidth=1,
            capsize=1, zorder=2)

# astericks and rsID (with dynamic upper bound of confidence interval)
texts = []
for _, r in df[df["sig"]!=""].iterrows():
    x  = r["POS"]
    y0 = r[upper_col]                # ← qui
    star_tr = offset_copy(ax.transData, fig, y=4, units='points')
    ax.text(x, y0, r["sig"], ha="center", va="bottom",
            transform=star_tr, color="firebrick",
            fontsize=7, zorder=4)
    rs_tr = offset_copy(ax.transData, fig, y=12, units='points')
    txt = ax.text(x, y0, r["rsID"], ha="center", va="bottom",
                  transform=rs_tr, color="black",
                  fontsize=7, zorder=4)
    texts.append(txt)
adjust_text(texts, arrowprops=dict(arrowstyle="-", color='grey', lw=.5))


ax.set_xlabel("Genomic position (bp)")
ax.set_ylabel("Selection coefficient (s)")
ax.set_title("* P<0.05   ** P<0.01   *** P<0.001")
ax.grid(alpha=.5)

ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x)}"))

cbar = fig.colorbar(sc, ax=ax, pad=0.01)
cbar.set_label("−log10 P  (0 – 10)")

fig.tight_layout()
fig.savefig(f"{prefix}_singleplot_ci.png", dpi=300)
fig.savefig(f"{prefix}_singleplot_ci.pdf")
plt.close(fig)
print("Plots written to:", Path(f"{prefix}_singleplot_ci.png/pdf").resolve())
PYCODE
  printf "\n"
  info "Phase‑2 completed."
echo
}
# ---------- PHASE-3  (re-inference + sweep-dating for ONE SNP) --------------
phase3() {
  PH="phase3"
  echo -e "\n          ╔═════════════════════════════════════════════════════════╗"
  echo -e "          ║        ⏱️  PHASE-3 – Dating a selective sweep  ⏱️         ║"
  echo -e "          ║   Please read the manual carefully before proceeding    ║"
  echo -e "          ╚═════════════════════════════════════════════════════════╝\n"
  # ------------------------ user input --------------------------------------
  read -rp "Choose chromosome to analyze (e.g. 2, 17, X): "  chr
  read -rp "Same prefix population name as used in Phase-2 (e.g. Finnish): "    pop
  read -rp "rsID of SNP to date   (e.g. rs123): "  rsid
  read -rp "df score for CLUES2 (e.g. 600): "      df_score

  # ------------------------ paths produced in phase-1/2 ---------------------
  P1_DIR="${WORK_BASE}/phase1/${pop}_chr${chr}"
  P2_DIR="${WORK_BASE}/phase2/${pop}_chr${chr}"
  [[ -d "$P1_DIR" && -d "$P2_DIR" ]] || \
        error "phase1/phase2 output not found for ${pop}_chr${chr}"
  
  # Phase-3 working directory
  RUN_ID="${pop}_chr${chr}"
  WORK_DIR="${WORK_BASE}/phase3/${RUN_ID}"
  mkdir -p "$WORK_DIR"
  printf "\n"
  info "[INFO]Phase-3 working dir: $WORK_DIR"

  COAL="${P1_DIR}/${pop}_EPS4COAL_chr${chr}.coal"
  [[ -f "$COAL" ]] || error "*.coal file missing"

  DERIVED="${P1_DIR}/${pop}_Derived_${rsid}.txt"
  [[ -f "$DERIVED" ]] || error "Derived file missing ($DERIVED)"

  # ----- derived-allele frequency (needed by inference.py) ------------------
  FREQ_FILE=$(ls "${P1_DIR}/${pop}_Frequency_chr${chr}_"*.txt | head -1)
  [[ -f "$FREQ_FILE" ]] || error "frequency file not found"
  freq=$(awk -v RSID="$rsid" '$1==RSID{print $4; exit}' "$FREQ_FILE")
  [[ -n "$freq" ]] || error "derived allele frequency for $rsid not found"

  # ------------------------ output directories ------------------------------
  DAT_DIR="${WORK_DIR}/Dating";              mkdir -p "$DAT_DIR"
  INF_OUT="${DAT_DIR}/${rsid}"             # prefix for first run
  OUT_JSON="${DAT_DIR}/${rsid}_Dating.json"

  RELATE_ROOT="$(dirname "$0")/Relate"
  CLUES_DIR="$(dirname "$0")/CLUES2"
  inf_py="${CLUES_DIR}/inference.py"
  rtc_py="${CLUES_DIR}/RelateToCLUES.py"

# -------------------------------------------------------------------------
# (A)  Sliding-window scan: one CLUES run per interval
# -------------------------------------------------------------------------
TIMES0="${P2_DIR}/${pop}_times_chr${chr}/${rsid}_times.txt"
[[ -f "$TIMES0" ]] || error "times file missing ($TIMES0)"

# Copiamo il file–times nella cartella di Phase-3 (così rimane tutto lì)
TIMES_FILE="${WORK_DIR}/${rsid}_times.txt"
cp "$TIMES0" "$TIMES_FILE"

# break-points:
BREAKS=(0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 900 1000 1200 1400 1600 1800 2000)

echo -e "\n▶  Scanning windows:"
WIN_PREFIXES=()

for (( i=0; i<${#BREAKS[@]}-1; i++ )); do
    left=${BREAKS[i]}
    right=${BREAKS[i+1]}
    out_pref="${INF_OUT}_${left}_${right}"
    WIN_PREFIXES+=("$out_pref")

    echo "   •  window [$left - $right]"

    if (( left == 0 )); then        # prima finestra: NIENTE --timeBins
        python "$inf_py" \
            --times   "$TIMES_FILE" \
            --popFreq "$freq" \
            --out     "$out_pref" \
            --df "$df_score" \
            --tCutoff "$right" \
            --coal    "$COAL" \
            --noAlleleTraj  >/dev/null
    else                            # finestre successive: un solo breakpoint
        python "$inf_py" \
            --times   "$TIMES_FILE" \
            --popFreq "$freq" \
            --out     "$out_pref" \
            --df "$df_score" \
            --tCutoff "$right" \
            --timeBins "$left" \
            --coal    "$COAL" \
            --noAlleleTraj  >/dev/null
    fi
done


# -------------------------------------------------------------------------
#  (A-2)  Raccogli i risultati dei singoli file e trova l’onset
# -------------------------------------------------------------------------
export INF_OUT rsid pop chr # let's read the variable for PY block
python - <<'PY'
import glob, os, re, pandas as pd, sys, json, pathlib, numpy as np

pref_root = os.environ['INF_OUT']         # es.  .../rs4988235
pattern   = f"{pref_root}_*_*.txt"        # match *_start_end_inference.txt
files     = sorted(glob.glob(pattern),
                   key=lambda p: int(re.search(r'_(\d+)_', p).group(1)),
                   reverse=True)          # ancient → recent

windows = []
for f in files:
    m = re.search(r'_(\d+)_(\d+)_inference\.txt$', f)
    if not m:  continue
    start,end = map(int, m.groups())
    row = pd.read_csv(f, sep="\t").iloc[0]
    # prendi SEMPRE l’ultima colonna SelectionMLE k
    last_k = max(int(c.split("SelectionMLE")[1])
                 for c in row.index if c.startswith("SelectionMLE"))
    s = row[f"SelectionMLE{last_k}"]
    p   = 10**(-row["-log10(p-value)"]) if "-log10(p-value)" in row else np.nan
    windows.append( (start,end,s,p) )

# ------------- criterio onset ---------------------------------------------
def first_trend(ws, k):
    for i in range(len(ws)-k+1):
        seg = ws[i:i+k]         # antico → recente
        ss  = [x[2] for x in seg]
        if all(val>0 for val in ss) and all(ss[j]>=ss[j-1] for j in range(1,k)):
            return seg[0]
    return None

onset = ( first_trend(windows,3) or
          first_trend(windows,2) or
          next((w for w in windows if w[2]>0), None) or
          max(windows, key=lambda t:t[2]) )

st,en,s,_ = onset
print(f"\nInitial onset ≈ {st} generations  (≈ {st*28} years)")
print(f"   epoch   : {st} – {en}")
print(f"   s(MLE)  : {s:.5f}\n")

# opzionale: salva un JSON intermedio
quick_out = pathlib.Path(f"{pref_root}_InitialOnset_Dating.json")
quick_out.write_text(json.dumps(dict(
        rsID=os.environ["rsid"],
        population        = os.environ["pop"],
        chromosome        = os.environ["chr"],
        onset_gen = int(st),
        onset_years = int(st*28),
        epoch_start = st,
        epoch_end   = en,
        epoch_start_years   = int(st*28),
        epoch_end_years = int(en*28),
        s_MLE       = round(s,5)
    ), indent=2))
PY

# -------------------------------------------------------------------------
#  (A-2)  initial onset done – JSON for first dating saved
# -------------------------------------------------------------------------

# >>> ask if the user want proceed with the bootstraps, otherwise exit <<<
read -rp "Proceed with bootstrap dating? [y/N]: " REPLY
if [[ ! "$REPLY" =~ ^[Yy]$ ]]; then
    echo -e "\nBootstrap step skipped – Phase-3 completed."
    return            # esce dalla funzione phase3 (o 'exit 0' se non è in funzione)
fi

# -----------------------------------------------------------------
#  (B)  bootstrap parameters
# -----------------------------------------------------------------
echo -e "\n▶  Bootstrap settings"
read -rp "Epochs to scan before initial onset (e.g. 200): "   G_START
read -rp "Epochs to scan after initial onset (e.g. 500): " G_END
read -rp "Non overlapping windows size  [default is 25]: "       STEP ; STEP=${STEP:-25}
read -rp "Number of bootstrap replicates  [default is 100]: "       NBOOT; NBOOT=${NBOOT:-100}
read -rp "df score for CLUES2  [default is 450]: "      DF   ; DF=${DF:-450}
read -rp "Importance sampling of branch lengths: " num_samples

# ---------- build BREAKS ----------------------------------------------------
if declare -F mapfile >/dev/null; then          # Bash ≥ 4
    mapfile -t BREAKS < <( seq "$G_START" "$STEP" "$G_END" )
else                                            # fallback portabile
    IFS=$'\n' read -r -d '' -a BREAKS < <( seq "$G_START" "$STEP" "$G_END" ; printf '\0' )
fi

last_break=${BREAKS[${#BREAKS[@]}-1]}
if (( last_break != G_END )); then
    BREAKS+=( "$G_END" )
fi

# quick sanity-check
(( ${#BREAKS[@]} == 0 )) && error "seq produced no break-points"

TCUT=$(( G_END + STEP ))        # --tCutoff passato a CLUES

echo -e "\nBreak-points : ${BREAKS[*]}"
echo   "tCutoff      : $TCUT"
echo   "Replicates   : $NBOOT"
echo

##############################################################################
#  C)  N BOOTSTRAP REPLICATES                           #
##############################################################################
BOOT_DIR="${DAT_DIR}/bootstrap_${rsid}"; mkdir -p "$BOOT_DIR"
GS_PREFIX="${P1_DIR}/${pop}_GS+COAL_chr${chr}"
POS=$(grep -m1 -w "$rsid" "${P1_DIR}/${pop}_SNPs_chr${chr}_"*.txt | cut -f2) \
     || error "bp position for $rsid not found"

echo -e "\n▶  Generating $NBOOT bootstrap trees & CLUES runs"

for (( rep=1; rep<=NBOOT; rep++ )); do
    seed=$(( ( RANDOM << 15 ) ^ RANDOM ))
    info "[bootstrap $rep]  SampleBranchLengths (seed=$seed)"
    nw_pref="${BOOT_DIR}/nw_${rep}"

    "${RELATE_ROOT}/scripts/SampleBranchLengths/SampleBranchLengths.sh" \
        --input "$GS_PREFIX" --output "$nw_pref" \
        --first_bp "$POS" --last_bp "$POS" \
        --format n --num_samples "$num_samples" --coal "$COAL" \
        --seed "$seed" -m 1.25e-8 &>/dev/null

    info "  RelateToCLUES"
    python "$rtc_py" --DerivedFile "$DERIVED" \
                     --RelateSamples "${nw_pref}.newick" \
                     --out "$nw_pref" &>/dev/null

# ---- CLUES runs each windows --------------------------------------
for (( i=0; i<${#BREAKS[@]}-1; i++ )); do
    left=${BREAKS[i]} ; right=${BREAKS[i+1]}
    win_pref="${nw_pref}_${left}_${right}"
    #log="${win_pref}_clues.log"
    
    info "    CLUES inference • window [$left - $right]"

    cmd=( python "$inf_py" --times "${nw_pref}_times.txt"
          --popFreq "$freq" --out "$win_pref"
          --df "$DF" --CI 0.95 --tCutoff "$right"
          --coal "$COAL" --noAlleleTraj )
    (( left > 0 )) && cmd+=( --timeBins "$left" )

    "${cmd[@]}" >/dev/null 2>&1 || {
        warn "      ➜ CLUES failed – window skipped"
        continue
    }
done

# ---------- AFTER all windows: concat Epoch-2 -------------------------
combo="${BOOT_DIR}/bootstrap_${rep}_${rsid}.txt"
: > "$combo"                            # svuota/crea

# header fisso
printf "logLR\t-log10(p-value)\tEpoch2_start\tEpoch2_end\tSelectionMLE2\n" \
       >> "$combo"

for f in "${nw_pref}"_*"_inference.txt"; do
    [[ -f $f ]] || continue

    awk -F'\t' '
        NR==1 {
            for(i=1;i<=NF;i++){
                if($i=="logLR")            c_l=i;
                if($i=="-log10(p-value)")  c_p=i;
                if($i=="Epoch2_start")     c_s=i;
                if($i=="Epoch2_end")       c_e=i;
                if($i=="SelectionMLE2")    c_m=i;
            }
            next
        }
        # stampa solo se tutte le colonne esistono nella riga
        (c_l && c_p && c_s && c_e && c_m && NF>=c_m) {
            printf("%s\t%s\t%s\t%s\t%s\n",
                   $c_l, $c_p, $c_s, $c_e, $c_m);
        }
    ' "$f" >> "$combo"
done

rm "${nw_pref}"_*"_inference.txt"



done   # -------------------------end bootstrap cycle----------------------


##############################################################################
#  D)  CALCULATE ONSET FOR ALL BOOTSTRAP                                    #
##############################################################################
output="onset_bootstraps.txt"
> "$output"   # svuota o crea il file

echo -e "\n▶  Computing onset for each bootstrap…"
onset_vals=()

for (( rep=1; rep<=NBOOT; rep++ )); do
    combo="${BOOT_DIR}/bootstrap_${rep}_${rsid}.txt"
    if [[ ! -s $combo ]]; then
        warn "Bootstrap $rep – file $combo empty or missing, skip"
        continue
    fi

    onset=$(python - <<'PY' "$combo"
import sys, pandas as pd

wins=[]
df=pd.read_csv(sys.argv[1], sep="\t")

for _, row in df.iterrows():
    st, en = int(row["Epoch2_start"]), int(row["Epoch2_end"])
    k  = max(int(c.split("SelectionMLE")[1]) for c in row.index if c.startswith("SelectionMLE"))
    s  = row[f"SelectionMLE{k}"]
    wins.append((st, en, s))

# ordine antico → recente
wins.sort(key=lambda w: w[0], reverse=True)

def streak(ws, k):
    for i in range(len(ws)-k+1):
        if all(w[2] > 0 for w in ws[i:i+k]):
            return ws[i][0]
    return None

o = streak(wins, 3) or streak(wins, 2) \
    or next((w[0] for w in wins if w[2] > 0), None) \
    or max(wins, key=lambda w: w[2])[0]

print(o)
PY
)

    # formato riga sia per terminale che per file
    line=$(printf "Bootstrap %2d : %s generations" "$rep" "$onset")
    echo  "$line"        # stampa a terminale
    echo  "$line" >> "$output"  # salva in append

    [[ $onset =~ ^[0-9]+$ ]] && onset_vals+=("$onset")
done



##############################################################################
#  E)  STATISTICHE                                                           #
##############################################################################
echo -e "\n▶  Summary statistics"
echo "Valid onsets collected: ${#onset_vals[@]}"

if (( ${#onset_vals[@]} < 2 )); then
    error "Fewer than 2 valid onsets – CI cannot be computed"
fi

stats=$(python3 - <<'PY' "${onset_vals[@]}"
import sys, numpy as np, math
vals = np.array(list(map(int, sys.argv[1:])), dtype=int)
print(int(np.median(vals)),
      int(math.floor(np.percentile(vals, 2.5))),
      int(math.ceil(np.percentile(vals, 97.5))))
PY
)

# --- split robusto, forziamo IFS a “spazio” / “newline” --------------------
IFS=$' \t\n' read -r MEDIAN CI_LOW CI_HIGH <<< "$stats"

echo
#echo "Median onset : $MEDIAN gen (~$(( MEDIAN * 28 )) yr)"
echo "95 % CI      : $CI_LOW – $CI_HIGH gen "\
     "(~$(( CI_LOW * 28 )) – $(( CI_HIGH * 28 )) yr)"

##############################################################################
#  F)  JSON OUTPUT                                                           #
##############################################################################
python3 - <<PY
import json, pathlib, os
out = dict(
    rsID              = os.environ["rsid"],
    population        = os.environ["pop"],
    chromosome        = os.environ["chr"],
    #onset_median_gen  = int("$MEDIAN"),
    #onset_median_year = int("$MEDIAN")*28,
    CI95_low_gen      = int("$CI_LOW"),
    CI95_high_gen     = int("$CI_HIGH"),
    CI95_low_year     = int("$CI_LOW")*28,
    CI95_high_year    = int("$CI_HIGH")*28,
    bootstraps        = int("$NBOOT"),
    step              = int("$STEP"),
    method            = "4 consecutive windows with s > 0"
)

out_dir = pathlib.Path("$BOOT_DIR").parent          # ← usa la variabile bash
out_path = out_dir / f"{os.environ['rsid']}_Boostraps_onset_Dating+CI.json"
out_path.write_text(json.dumps(out, indent=2))

print(f"JSON saved → {out_path}")
PY

}
# ---------- dispatch ----------------------------------------------------------
case "$OPTION" in
  1) phase1 ;;
  2) phase2 ;;
  3) phase3 ;;
  *) error "Invalid option" ;;
esac

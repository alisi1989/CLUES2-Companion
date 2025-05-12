#!/usr/bin/env python3
# plot_traj.py  – visualizza la heat‑map delle posterior frequency trajectories

import argparse
import gzip
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colormaps   # <‑ nuovo import (Matplotlib ≥ 3.7)

# ----------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--freqs',  required=True, type=str,
                    help='CSV con la griglia di frequenze (una colonna).')
parser.add_argument('--post',   required=True, type=str,
                    help='CSV con la matrice di posterior log‑prob (freq x epoch).')
parser.add_argument('--figure', required=True, type=str,
                    help='Prefisso del file di figura in output (senza estensione).')
parser.add_argument('--posterior_intervals', default=[0.5, 0.75, 0.95, 0.999],
                    type=float, nargs='+',
                    help='Posterior thresholds per i diversi intervalli di confidenza.')
parser.add_argument('--generation_time', default=-1.0, type=float,
                    help='Anni per generazione (se >0 converte l’asse x in anni).')
args = parser.parse_args()

# ----------------------------------------------------------------------
# prepara intervalli e colori
ConfInts = sorted(round(x, 4) for x in args.posterior_intervals)
ColorInt = [1.0 - i / len(ConfInts) for i in range(len(ConfInts))]

freqs    = np.loadtxt(args.freqs, delimiter=",", dtype=float)
logpost  = np.loadtxt(args.post, delimiter=",", dtype=float)
epochs   = np.linspace(0, logpost.shape[1], logpost.shape[1] + 1)

# matrice da plottare
MATRIX   = np.zeros_like(logpost)

# --- trova indice di partenza ragionevole per ogni colonna
StartIdx = [0] * (len(epochs) - 1)
for t in range(len(epochs) - 1):
    slice_t = logpost[:, t]
    nz      = np.where(slice_t > 1e-6 / len(freqs))[0]
    if nz.size:
        StartIdx[t] = nz[0]

# --- costruisci la heat‑map per ciascun intervallo di confidenza
for k, ci in enumerate(ConfInts):
    for t in range(len(epochs) - 1):
        best_span = np.inf
        lo = hi = -1
        slice_t = logpost[:, t]
        for i in range(StartIdx[t], len(freqs)):
            p_sum = slice_t[i]
            for j in range(i + 1, len(freqs)):
                p_sum += slice_t[j]
                span = freqs[j] - freqs[i]
                if p_sum >= ci and span < best_span:
                    best_span, lo, hi = span, i, j
                elif span > best_span:
                    break
        if lo >= 0:
            MATRIX[lo:hi + 1, t] = np.maximum(MATRIX[lo:hi + 1, t], ColorInt[k])

# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(20, 10))
x_vals  = epochs[:-1] if args.generation_time < 0 else epochs[:-1] * args.generation_time
x_label = ('Generations before present' if args.generation_time < 0
           else 'Years before present')

pcm = ax.pcolormesh(x_vals, freqs, MATRIX, shading='auto')
ax.set_xlim(x_vals[0], x_vals[-1])
ax.set_ylim(0, 1.0)
ax.set_xlabel(x_label, fontsize=20)
ax.set_ylabel('Derived allele frequency', fontsize=20)
ax.tick_params(labelsize=18)

# legenda
cmap     = colormaps.get_cmap('viridis')           # <‑ nuovo modo
handles  = [mpatches.Patch(color=cmap(ci), label=f'{ci*100:.1f}% posterior interval')
            for ci in ColorInt]
ax.legend(handles=handles, loc='upper right')

# salvataggio figura
plt.savefig(f'{args.figure}.png', format='png', dpi=300, bbox_inches='tight')
plt.close()

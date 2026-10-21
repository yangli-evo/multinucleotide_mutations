# Multinucleotide mutations

Code accompanying Li & Zhang, *Multinucleotide mutations expand evolutionary accessibility on rugged adaptive landscapes*.

About 3% of point mutations are multinucleotide mutations (MNMs), which change two or more nearby nucleotides at the same time. On rugged fitness landscapes, MNMs let adaptive walks escape local fitness peaks. This repository contains a simulator of adaptive walks with single-nucleotide mutations (SNMs) and MNMs, and the scripts used for the analyses in the paper.

## Features

- Generate single-nucleotide mutants (SNMs), adjacent double mutants (adjacent MNMs; sites *i* and *i*+1) and all double mutants (general MNMs; any two sites *i* < *j*).
- Represent complete or incomplete fitness landscapes over any alphabet (DNA, amino acids, two-state sites). Genotypes missing from a landscape are treated as non-viable.
- Simulate adaptive walks with SNMs only, or with SNMs plus adjacent or general MNMs (3% MNMs by default), under
  - the **strong-selection, weak-mutation (SSWM)** scheme, where a random beneficial mutation is accepted, and
  - the **greedy** scheme, where the fittest beneficial mutation is accepted.
- Command-line program `mnm-simulate` that runs on any landscape given as a table.

The walk rule follows the Methods of the paper. If beneficial SNMs and beneficial MNMs both exist, the type of the next mutation is drawn according to the mutation proportions (97% SNMs, 3% MNMs), and a beneficial mutation of that type is accepted. At a local peak (no beneficial SNM), a beneficial MNM is accepted if one exists. The walk ends at a genotype from which no mutation of the allowed classes increases fitness.

## Install

```bash
pip install -e .                 # core package (numpy, pandas)
pip install -e ".[analysis]"     # + scipy, matplotlib and jupyter for the analysis scripts and notebook
```

## Python

```python
from mnm_mutations import (
    generate_single_base_mutants,
    generate_adjacent_double_mutants,
    generate_double_mutants,
)

seq = "ATCG"

snm = generate_single_base_mutants(seq)
adj_mnm = generate_adjacent_double_mutants(seq)
mnm = generate_double_mutants(seq)

print(len(snm))     # 12 = 3L
print(len(adj_mnm)) # 27 = 9(L-1)
print(len(mnm))     # 54 = 9L(L-1)/2
```

Adaptive walks on a landscape:

```python
import numpy as np
from mnm_mutations import load_table, walk_original, peak_thresholds

ls = load_table("examples/hoc_L5_example.txt")        # columns: seq_info, fitness
rng = np.random.RandomState(1)
start = ls.encode("ACGTA")

path_snm = walk_original(ls, start, None, rng=rng)              # SNMs only
path_adj = walk_original(ls, start, "adj", p=0.03, rng=rng)     # + 3% adjacent MNMs
path_gen = walk_original(ls, start, "gen", p=0.03, rng=rng)     # + 3% general MNMs
path_greedy = walk_original(ls, start, "gen", p=0.03, rng=rng, greedy=True)
print([ls.decode(g) for g in path_gen], ls.F[path_gen[-1]])

# fitness of the global peak and threshold for high peaks (top 10% of peaks)
global_fitness, high_threshold = peak_thresholds(ls.F[ls.snm_peaks()])
```

## Command line

```bash
# SSWM walks with 3% adjacent MNMs from 3,000 random starting genotypes
mnm-simulate landscape.txt --mnm adj --n-starts 3000 --seed 1 -o adj.tsv
# SNM-only and general-MNM walks under the greedy scheme
mnm-simulate landscape.txt --mnm none --greedy -o snm_greedy.tsv
mnm-simulate landscape.txt --mnm gen  --greedy -o gen_greedy.tsv
# 9% instead of 3% MNMs
mnm-simulate landscape.txt --mnm gen --mnm-fraction 0.09 -o gen_9pct.tsv
# every genotype as a start, all four nucleotides possible at every site (TF-binding landscapes)
mnm-simulate tf.txt --mnm adj --all-starts --alphabet ACGT -o tf_adj.tsv
```

Input: a tab-separated table with columns `seq_info` and `fitness`. Output: one row per genotype visited (`mutant_step`, `seq_info`, `fitness`, `start_sequence`). Run `mnm-simulate --help` for all options.

## Reproducing the analyses

The landscapes and simulated trajectories are archived at Zenodo (see [`data/README.md`](data/README.md)). Set `MNM_DATA` to the unpacked archive. The analysis scripts need Python ≥ 3.9 with scipy; the figure scripts need R ≥ 4.3 with ggplot2, ggpubr, dplyr, tidyr and patchwork.


```bash
cd analysis
python 01_mnm_proportion_walks.py      # walks with 3%, 9% and 27% MNMs on all landscapes (~15 min on 14 cores)
python 02_escape_vs_exploration.py     # exploration (y - x) and escape (z) per landscape and landscape class
python 03_peak_cv.py                   # CV of peak height and peak probabilities of the TF-binding landscapes
Rscript R/fig_S19_peak_cv.R
Rscript R/fig_S20_escape_exploration.R
```

The tables plotted in Figs. S19 and S20 are included in [`results/`](results/), so the two R scripts run without the Zenodo data. Figures are written to `figures/`.

## Repository layout

```
src/mnm_mutations/   mutants.py (mutant generation), landscape.py (landscapes, adaptive walks), simulate.py (CLI)
tests/               pytest suite
examples/            a small House-of-Cards landscape (L = 5) for testing
notebooks/           Fig1_HoC.ipynb (Fig. 1 and fitness-distribution dependence of the HoC results)
original_pipeline/   scripts that simulated the adaptive walks and produced Figs. 1-4, S1-S18, S21-S25
analysis/            scripts for Figs. S19 and S20 (Python) and their figures (analysis/R)
results/             tables produced by the analysis scripts
data/                placeholder for the Zenodo data (see data/README.md)
```

## Tests

```bash
pip install pytest
pytest -q
```

## Citation

Li, Y. & Zhang, J. Multinucleotide mutations expand evolutionary accessibility on rugged adaptive landscapes. *Nature Communications* (2026).

## License

MIT (see `LICENSE`).

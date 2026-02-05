# multinucleotide_mutations

Exploring the role of multi-mutation leaps to reach a global peak.

In many datasets, multinucleotide mutations (MNMs) account for only ~3% of all nucleotide substitutions, but they can disproportionately affect accessibility of high-fitness genotypes.

## Features

This repo provides simple utilities to generate:
- Single-nucleotide mutants (SNM)
- Adjacent double mutants (adjacent MNM; positions i and i+1)
- All double mutants (any pair of positions i < j)

## Install

```bash
pip install -e .

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

print(len(snm))     # 12 = 4*3
print(len(adj_mnm)) # 27 = (4-1)*9
print(len(mnm))     # 54 = C(4,2)*9

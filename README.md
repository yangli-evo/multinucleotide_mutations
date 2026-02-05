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


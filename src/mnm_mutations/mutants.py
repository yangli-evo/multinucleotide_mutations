# src/mnm_mutations/mutants.py
from __future__ import annotations

from itertools import combinations, product
from typing import List, Sequence


_BASES: Sequence[str] = ("A", "T", "C", "G")


def validate_dna_sequence(dna_sequence: str) -> str:
    """
    Validate and normalize a DNA sequence.

    Rules:
    - Must be a non-empty string
    - Allowed characters: A/T/C/G (case-insensitive)
    """
    if not isinstance(dna_sequence, str):
        raise TypeError(f"dna_sequence must be a str, got {type(dna_sequence)}")

    seq = dna_sequence.strip().upper()
    if len(seq) == 0:
        raise ValueError("dna_sequence must be non-empty after stripping whitespace")

    bad = sorted(set([c for c in seq if c not in _BASES]))
    if bad:
        raise ValueError(
            f"Invalid character(s) in dna_sequence: {bad}. Allowed: A/T/C/G only."
        )
    return seq


def generate_single_base_mutants(dna_sequence: str, *, deduplicate: bool = False) -> List[str]:
    """
    Generate all single-nucleotide substitution mutants for a DNA sequence.

    For a sequence of length n, returns n * 3 mutants.
    """
    seq = validate_dna_sequence(dna_sequence)
    mutants: List[str] = []

    for i, current_base in enumerate(seq):
        for base in _BASES:
            if base != current_base:
                mutant = seq[:i] + base + seq[i + 1 :]
                mutants.append(mutant)

    if deduplicate:
        # preserve order
        seen = set()
        mutants = [m for m in mutants if not (m in seen or seen.add(m))]
    return mutants


def generate_adjacent_double_mutants(dna_sequence: str, *, deduplicate: bool = False) -> List[str]:
    """
    Generate all adjacent double-substitution mutants (i, i+1).

    For a sequence of length n, returns (n-1) * 3 * 3 mutants.
    """
    seq = validate_dna_sequence(dna_sequence)
    mutants: List[str] = []

    for i in range(len(seq) - 1):
        b1 = seq[i]
        b2 = seq[i + 1]
        for new_b1 in _BASES:
            if new_b1 == b1:
                continue
            for new_b2 in _BASES:
                if new_b2 == b2:
                    continue
                mutant = seq[:i] + new_b1 + new_b2 + seq[i + 2 :]
                mutants.append(mutant)

    if deduplicate:
        seen = set()
        mutants = [m for m in mutants if not (m in seen or seen.add(m))]
    return mutants


def generate_double_mutants(dna_sequence: str, *, deduplicate: bool = False) -> List[str]:
    """
    Generate all (non-adjacent + adjacent) double-substitution mutants across all pairs (i, j), i<j.

    For a sequence of length n, returns C(n,2) * 3 * 3 mutants.
    """
    seq = validate_dna_sequence(dna_sequence)
    mutants: List[str] = []
    seq_list = list(seq)
    n = len(seq_list)

    for i, j in combinations(range(n), 2):
        options_i = [b for b in _BASES if b != seq_list[i]]
        options_j = [b for b in _BASES if b != seq_list[j]]
        for bi, bj in product(options_i, options_j):
            new_seq = seq_list.copy()
            new_seq[i] = bi
            new_seq[j] = bj
            mutants.append("".join(new_seq))

    if deduplicate:
        seen = set()
        mutants = [m for m in mutants if not (m in seen or seen.add(m))]
    return mutants

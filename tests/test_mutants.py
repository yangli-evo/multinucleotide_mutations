# tests/test_mutants.py
import pytest

from mnm_mutations import (
    generate_adjacent_double_mutants,
    generate_double_mutants,
    generate_single_base_mutants,
    validate_dna_sequence,
)


def test_validate_dna_sequence_ok():
    assert validate_dna_sequence("atcg") == "ATCG"


def test_validate_dna_sequence_bad_char():
    with pytest.raises(ValueError):
        validate_dna_sequence("ATNG")


def test_single_mutants_count():
    # length 4 -> 4*3 = 12
    muts = generate_single_base_mutants("ATCG")
    assert len(muts) == 12
    assert "TTCG" in muts  # A->T at pos0
    assert "ATCA" in muts  # G->A at pos3


def test_adjacent_double_mutants_count():
    # length 4 -> (4-1)*9 = 27
    muts = generate_adjacent_double_mutants("ATCG")
    assert len(muts) == 27


def test_double_mutants_count():
    # length 4 -> C(4,2)*9 = 6*9 = 54
    muts = generate_double_mutants("ATCG")
    assert len(muts) == 54


def test_deduplicate_preserves_size_for_normal_cases():
    # For standard inputs there should be no duplicates anyway, but keep the behavior stable
    muts1 = generate_double_mutants("ATCG", deduplicate=False)
    muts2 = generate_double_mutants("ATCG", deduplicate=True)
    assert len(muts1) == len(muts2)
    assert muts1 == muts2

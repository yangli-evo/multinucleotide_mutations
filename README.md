# multinucleotide_mutations
Exploring the role of multimutation leaps to reach global peak multinucleotide mutations (MNMs) account for only ~3% of all nucleotide substitutions

## Install

```bash
pip install -e .

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


---

## 8) add this command to your repo

```bash
# under the repo
mkdir -p src/mnm_mutations tests .github/workflows

# copy and save below python file：
# src/mnm_mutations/mutants.py
# src/mnm_mutations/__init__.py
# tests/test_mutants.py
# pyproject.toml
# .github/workflows/tests.yml

pip install -e .
pytest -q

git add .
git commit -m "Add mutant generation utilities (SNM/adjacent MNM/double MNM) with tests and CI"
git push

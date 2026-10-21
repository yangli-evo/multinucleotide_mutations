# src/mnm_mutations/simulate.py
import argparse
import numpy as np
import pandas as pd
from .landscape import Landscape, walk_original


def main(argv=None):
    ap = argparse.ArgumentParser(prog='mnm-simulate', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('landscape')
    ap.add_argument('--mnm', choices=['none', 'adj', 'gen'], default='none')
    ap.add_argument('--mnm-fraction', type=float, default=0.03)
    ap.add_argument('--greedy', action='store_true')
    ap.add_argument('--n-starts', type=int, default=3000)
    ap.add_argument('--all-starts', action='store_true', help='start from every genotype in the landscape')
    ap.add_argument('--starts-file', help='file with one starting genotype per line')
    ap.add_argument('--alphabet', help="states per site, e.g. ACGT (default: inferred per site)")
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('-o', '--out', required=True)
    a = ap.parse_args(argv)

    rng = np.random.RandomState(a.seed)
    df = pd.read_table(a.landscape)
    ls = Landscape(df.seq_info.astype(str).values, df.fitness.values, a.landscape, a.alphabet)
    if a.starts_file:
        starts = [ls.encode(s.strip()) for s in open(a.starts_file) if s.strip()]
    elif a.all_starts or a.n_starts >= len(ls.codes):
        starts = list(ls.codes)
    else:   # uniform sampling of starting genotypes without replacement
        starts = list(rng.choice(ls.codes, a.n_starts, replace=False))
    kind = None if a.mnm == 'none' else a.mnm
    rows = []
    for s in starts:
        path = walk_original(ls, s, kind, a.mnm_fraction, rng, greedy=a.greedy)
        st = ls.decode(s)
        rows += [(i, ls.decode(c), ls.F[c], st) for i, c in enumerate(path)]
    out = pd.DataFrame(rows, columns=['mutant_step', 'seq_info', 'fitness', 'start_sequence'])
    out.to_csv(a.out, sep='\t', index=False)


if __name__ == '__main__':
    main()

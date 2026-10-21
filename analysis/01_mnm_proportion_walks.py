import os
import sys
import numpy as np
import pandas as pd
from multiprocessing import Pool
import mnmlib as M

PS = [0.03, 0.09, 0.27]
OUT = os.path.join(M.RES, 'mnm_proportion')


def starts_of(group, tag):
    t = pd.read_table(M.traj_path(group, tag, 'snm'), usecols=['start_sequence'])
    return t.start_sequence.unique()


def run(args):
    group, tag, seed = args
    ls = M.load_landscape(group, tag)
    peaks = ls.snm_peaks()
    is_peak = np.zeros(ls.size, dtype=bool)
    is_peak[peaks] = True
    gmax, thr = M.peak_thresholds(ls.F[peaks])
    fs = np.sort(ls.F[ls.codes])
    q = lambda f: np.searchsorted(fs, f, 'right') / len(fs)
    W = np.array(ls.w)
    K = np.array(ls.k)
    starts = [ls.encode(s) for s in starts_of(group, tag)]
    rows = []
    for j, (kind, p) in enumerate([(None, 0.0)] + [(k, p) for k in ('adj', 'gen') for p in PS]):
        rng = np.random.RandomState(seed * 100 + j)
        for s in starts:
            path = np.array(M.walk_original(ls, s, kind, p, rng), dtype=np.int64)
            f = ls.F[path]
            pk = is_peak[path]
            dig = (path[:, None] // W[None, :]) % K[None, :]
            ham = (dig[1:] != dig[:-1]).sum(1)
            esc = pk[:-1]                                   # step leaves an SNM peak
            i1 = int(np.argmax(pk)) if pk.any() else len(path) - 1
            n = len(path) - 1
            sq, eq = q(f[0]), q(f[-1])
            rows.append((group, tag, kind or 'snm', p, n, int(esc.sum()), int(((ham == 2) & ~esc).sum()),
                         f[0], f[-1], sq, eq, (eq - sq) / n if n else np.nan,
                         (f[-1] - f[0]) / n if n else np.nan, q(f[i1]), int(pk.sum()),
                         bool(f[-1] >= gmax), bool(f[-1] >= thr)))
    cols = ['group', 'landscape', 'cls', 'p', 'n_steps', 'n_esc', 'n_expl', 'start_fit', 'end_fit', 'start_q', 'end_q',
            'step_q', 'step_fit', 'trap1_q', 'n_peaks_vis', 'global', 'high']
    df = pd.DataFrame(rows, columns=cols)
    esc_frac = {k: np.mean([len(ls.beneficial(int(c), k)) > 0 for c in peaks]) for k in ('adj', 'gen')}
    meta = dict(group=group, landscape=tag, n_genotypes=len(fs), n_peaks=len(peaks),
                peak_fraction=len(peaks) / len(fs), esc_frac_adj=esc_frac['adj'], esc_frac_gen=esc_frac['gen'])
    return df, meta


def summarise(df, meta):
    g = df.groupby(['group', 'landscape', 'cls', 'p'])
    s = g.agg(n_walks=('global', 'size'), p_global=('global', 'mean'), p_high=('high', 'mean'),
              mean_end_q=('end_q', 'mean'), mean_end_fit=('end_fit', 'mean'), mean_steps=('n_steps', 'mean'),
              mean_step_q=('step_q', 'mean'), mean_step_fit=('step_fit', 'mean'), mean_trap1_q=('trap1_q', 'mean'),
              mean_esc=('n_esc', 'mean'), mean_expl=('n_expl', 'mean'),
              frac_walks_mnm=('n_esc', lambda x: np.nan), mean_peaks_vis=('n_peaks_vis', 'mean')).reset_index()
    a = df.assign(any_mnm=(df.n_esc + df.n_expl) > 0, any_esc=df.n_esc > 0, any_expl=df.n_expl > 0)
    f = a.groupby(['group', 'landscape', 'cls', 'p'])[['any_mnm', 'any_esc', 'any_expl']].mean().reset_index()
    s = s.drop(columns='frac_walks_mnm').merge(f, on=['group', 'landscape', 'cls', 'p'])
    # escape vs exploration: benefit relative to the SNM-only walks from the same starts
    snm = s[s.cls == 'snm'][['group', 'landscape', 'mean_end_q']].rename(columns={'mean_end_q': 'snm_end_q'})
    s = s.merge(snm, on=['group', 'landscape'])
    s['exploration'] = s.mean_trap1_q - s.snm_end_q
    s['escape'] = s.mean_end_q - s.mean_trap1_q
    # fraction of the SNM peaks reached by the walks that were escaped (= escapable); 1 - this = stopped at a peak
    return s.merge(pd.DataFrame(meta), on=['group', 'landscape'])


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    groups = sys.argv[1:] or ['HoC', 'RMF', 'NK', 'MBE', 'TF']
    summ = []
    for group in groups:
        jobs = [(group, t, 1000 * (i + 1) + len(group)) for i, t in enumerate(M.tags(group))]
        with Pool(min(14, len(jobs))) as pool:
            res = pool.map(run, jobs, chunksize=1 if group != 'TF' else 8)
        df = pd.concat([r[0] for r in res], ignore_index=True)
        meta = [r[1] for r in res]
        df.to_csv(os.path.join(OUT, f'walks_{group}.tsv.gz'), sep='\t', index=False)
        s = summarise(df, meta)
        s.to_csv(os.path.join(OUT, f'summary_{group}.tsv'), sep='\t', index=False)
        print(group, 'done', len(df), flush=True)
        print(s[s.landscape == s.landscape.iloc[0]][['cls', 'p', 'p_high', 'p_global', 'mean_end_q', 'mean_steps',
                                                      'mean_esc', 'mean_expl', 'exploration', 'escape']].round(4).to_string(),
              flush=True)

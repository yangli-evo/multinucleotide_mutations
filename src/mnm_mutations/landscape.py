# src/mnm_mutations/landscape.py
import numpy as np
import pandas as pd

__all__ = ['Landscape', 'load_table', 'walk_original', 'peak_thresholds']


class Landscape:
    def __init__(self, seqs, fitness, name='', alphabet=None):
        seqs = list(seqs)
        self.name = name
        self.L = len(seqs[0])
        self.alpha = ([sorted(alphabet)] * self.L if alphabet else
                      [sorted(set(s[i] for s in seqs)) for i in range(self.L)])
        self.k = np.array([len(a) for a in self.alpha])
        self.w = np.ones(self.L, dtype=np.int64)
        for i in range(self.L - 2, -1, -1):
            self.w[i] = self.w[i + 1] * self.k[i + 1]
        self.size = int(np.prod(self.k))
        self.idx = [{c: j for j, c in enumerate(a)} for a in self.alpha]
        codes = np.array([self.encode(s) for s in seqs], dtype=np.int64)
        self.F = np.full(self.size, np.nan)
        self.F[codes] = np.asarray(fitness, dtype=float)
        self.codes = np.sort(codes)
        # shifts[i][d] = code shifts that move site i from state d to every other state
        self.shifts = [[np.array([(n - d) * self.w[i] for n in range(self.k[i]) if n != d], dtype=np.int64)
                        for d in range(self.k[i])] for i in range(self.L)]
        km1 = self.k - 1
        self.n_snm = int(km1.sum())
        self.n_adj = int(sum(km1[i] * km1[i + 1] for i in range(self.L - 1)))
        self.n_gen = int(sum(km1[i] * km1[j] for i in range(self.L) for j in range(i + 1, self.L)))
        self.pairs = {'adj': [(i, i + 1) for i in range(self.L - 1)],
                      'gen': [(i, j) for i in range(self.L) for j in range(i + 1, self.L)]}
        self._cache = {}

    # ---------- encoding ----------
    def encode(self, s):
        return int(sum(self.idx[i][c] * self.w[i] for i, c in enumerate(s)))

    def decode(self, code):
        d = self.digits(code)
        return ''.join(self.alpha[i][d[i]] for i in range(self.L))

    def digits(self, code):
        return [(code // self.w[i]) % self.k[i] for i in range(self.L)]

    # ---------- neighbours ----------
    def neighbors(self, code, kind):
        """Codes of all mutants of class `kind` ('snm', 'adj' or 'gen'), including absent genotypes."""
        key = (code, kind)
        if key in self._cache:
            return self._cache[key]
        d = self.digits(code)
        if kind == 'snm':
            out = code + np.concatenate([self.shifts[i][d[i]] for i in range(self.L)])
        else:
            out = np.concatenate([(self.shifts[i][d[i]][:, None] + self.shifts[j][d[j]][None, :]).ravel()
                                  for i, j in self.pairs[kind]]) + code
        if len(self._cache) > 2_000_000:
            self._cache.clear()
        self._cache[key] = out
        return out

    def beneficial(self, code, kind):
        """Codes of the mutants of class `kind` that are strictly fitter than `code`."""
        nb = self.neighbors(code, kind)
        f = self.F[nb]
        with np.errstate(invalid='ignore'):
            return nb[f > self.F[code]]

    def total(self, kind):
        """Total number of possible mutations of class `kind` (present or absent in the landscape)."""
        return {'snm': self.n_snm, 'adj': self.n_adj, 'gen': self.n_gen}[kind]

    # ---------- landscape-wide properties ----------
    def snm_peaks(self):
        """Codes of the peaks: genotypes strictly fitter than all existing single-mutant neighbours."""
        c = self.codes
        f = self.F[c]
        is_peak = np.ones(len(c), dtype=bool)
        dig = np.stack([(c // self.w[i]) % self.k[i] for i in range(self.L)], axis=1)
        for i in range(self.L):
            for n in range(self.k[i]):
                nb = c + (n - dig[:, i]) * self.w[i]
                ok = dig[:, i] != n
                fn = self.F[nb[ok]]
                with np.errstate(invalid='ignore'):
                    bad = fn >= f[ok]
                tmp = is_peak[ok]
                tmp[bad] = False
                is_peak[ok] = tmp
        return c[is_peak]


def load_table(path, name='', alphabet=None):
    """Read a tab-separated landscape with columns `seq_info` (genotype) and `fitness`."""
    df = pd.read_table(path, sep='\t')
    return Landscape(df['seq_info'].astype(str).values, df['fitness'].values, name, alphabet)


# ---------------- adaptive walks ----------------
def walk_original(ls, start, kind, p=0.03, rng=np.random, greedy=False):
    pick = (lambda c: c[np.argmax(ls.F[c])]) if greedy else (lambda c: c[rng.randint(len(c))])
    path = [start]
    g = start
    while True:
        bs = ls.beneficial(g, 'snm')
        if kind is None:
            if len(bs) == 0:
                break
            g = pick(bs)
        else:
            bm = ls.beneficial(g, kind)
            if len(bs) > 0:
                if len(bm) > 0 and rng.random_sample() >= 1 - p:
                    g = pick(bm)
                else:
                    g = pick(bs)
            elif len(bm) > 0:
                g = pick(bm)
            else:
                break
        path.append(g)
    return path


def peak_thresholds(peak_fitness):
    """Fitness of the global peak and the threshold defining high peaks (top 10% of peaks)."""
    pf = pd.Series(peak_fitness)
    return pf.max(), pf.quantile(0.90)

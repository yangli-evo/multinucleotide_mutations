import glob
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
try:
    from mnm_mutations.landscape import Landscape, load_table, walk_original, peak_thresholds
except ImportError:   # source checkout without `pip install -e .`
    sys.path.insert(0, os.path.join(ROOT, 'src'))
    from mnm_mutations.landscape import Landscape, load_table, walk_original, peak_thresholds

DATA = os.environ.get('MNM_DATA', os.path.join(ROOT, 'data'))
RES = os.path.join(ROOT, 'results')
FIG = os.path.join(ROOT, 'figures')

HOC_L = [5, 6, 7, 8, 9, 10]
RMF_SIGMA = ['0.5', '0.53', '0.56', '0.59', '0.62', '0.67', '0.71', '0.77', '0.83', '0.91', '1.0', '1.11',
             '1.25', '1.43', '1.67', '2.0', '2.5', '3.33', '5.0', '10.0']
NK_K = [1, 2, 3, 4, 5, 6, 7, 8]
MBE = ['E.coli_DHFR', 'E.coli_Shine-Dalgamo', 'tRNA', 'ACE2_affinity', 'Fluorescent_protein', 'GB1_protein']
MBE_LABEL = {'E.coli_DHFR': 'DHFR', 'E.coli_Shine-Dalgamo': 'Shine-Dalgarno', 'tRNA': 'tRNA',
             'ACE2_affinity': 'ACE2 affinity', 'Fluorescent_protein': 'Fluorescent protein',
             'GB1_protein': 'GB1'}

# simulated trajectories: walk class -> (sub-folder, file suffix)
TRAJ = {
    'snm': ('SNM', '_adaptive_walk_by_SNM_allsites_3k_trajectory.txt'),
    'adj': ('MNM', '_adaptive_walk_by_MNM_adjacent_allsites_3k_trajectory.txt'),
    'gen': ('MNM2site', '_adaptive_walk_by_MNM2sites_allsites_3k_trajectory.txt'),
}


def landscape_path(group, tag):
    if group == 'HoC':
        return f'{DATA}/01_HouseCard/{tag}.txt'
    if group == 'RMF':
        return f'{DATA}/02_RMF/SNM/{tag}.txt'
    if group == 'NK':
        return f'{DATA}/03_NK/01_updated_NK/000_landscape_peak/{tag}.txt'
    if group == 'MBE':
        return f'{DATA}/05_MBE/{tag}.txt'
    raise ValueError(group)


def group_dir(group, greedy=False):
    if greedy:   # greedy walks were simulated on the HoC landscapes
        return {'HoC': f'{DATA}/01_HouseCard_greedy'}[group]
    d = {'HoC': '01_HouseCard', 'RMF': '02_RMF', 'NK': '03_NK/01_updated_NK', 'TF': '04_TF', 'MBE': '05_MBE'}[group]
    return f'{DATA}/{d}'


def tags(group):
    if group == 'HoC':
        return [f'seq_{L}bp' for L in HOC_L]
    if group == 'RMF':
        return [f'RMF_9bp_std{s}' for s in RMF_SIGMA]
    if group == 'NK':
        return [f'NK_model_9bp_K{k}' for k in NK_K]
    if group == 'MBE':
        return MBE
    if group == 'TF':
        fs = glob.glob(f'{DATA}/04_TF/SNM/*' + TRAJ['snm'][1])
        return sorted(os.path.basename(f).replace(TRAJ['snm'][1], '') for f in fs)
    raise ValueError(group)


def traj_path(group, tag, cls, greedy=False):
    sub, suf = TRAJ[cls]
    if greedy and cls != 'snm':
        suf = suf.replace('_3k_trajectory', '_3k_greedy_trajectory')
    return f'{group_dir(group, greedy)}/{sub}/{tag}{suf}'


def load_landscape(group, tag):
    if group == 'TF':
        t = pd.read_table(traj_path('TF', tag, 'snm'), usecols=['mutant_step', 'seq_info', 'fitness'])
        t = t[t.mutant_step == 0].drop_duplicates('seq_info')
        return Landscape(t.seq_info.values, t.fitness.values, tag, alphabet='ACGT')
    return load_table(landscape_path(group, tag), tag)


def read_traj(path):
    t = pd.read_table(path, usecols=['mutant_step', 'seq_info', 'fitness', 'start_sequence'])
    return t.sort_values(['start_sequence', 'mutant_step'], kind='mergesort')

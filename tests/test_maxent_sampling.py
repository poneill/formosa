from maxent_sampling import *
from formosa_utils import mean_ci

def test_maxent_motif():
    trials = 1000
    N, L, desired_ic = 10,10,10
    motifs = [maxent_motif(N,L,desired_ic) for i in trange(trials)]
    lb, ub = mean_ci(map(motif_ic,motifs))
    assert lb < desired_ic < ub

def test_spoof_maxent_motifs():
    trials = 1000
    motif = ['CGGTGAACTA',
             'CGGTGTGCGA',
             'CGCTGTGCTG',
             'CGGGATGCAA',
             'CACGCTACGA',
             'CGCTATGCTA',
             'CGGTTGGCTA',
             'CGGCGTGCTA',
             'CGGTATATTG',
             'CGGGTTGCGA']
    given_ic = motif_ic(motif) # ~ 9.05 bits
    motifs = spoof_maxent_motifs(motif,trials)
    lb, ub = mean_ci(map(motif_ic,motifs))
    assert lb < given_ic < ub

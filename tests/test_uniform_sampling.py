from uniform_sampling import *

def test_uniform_motif():
    motif = uniform_motif(N=10,L=10,desired_ic=10,epsilon=0.1)
    assert 9.9 < motif_ic(motif) < 10.1

def test_spoof_motifs_uniform():
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
    spoofs = spoof_uniform_motifs(motif,trials)
    assert all(given_ic - 0.1 < motif_ic(spoof) < given_ic + 0.1 for spoof in spoofs)

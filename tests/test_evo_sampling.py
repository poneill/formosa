from formosa_utils import motif_ic, mean_ci, occupancies
from evo_sampling import *
from nose.tools import assert_less_equal

def test_spoof_motif_cftp():
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
    motifs = spoof_motif_cftp(motif, num_motifs=1000)
    spoof_ics = map(motif_ic, motifs)
    lb, ub = mean_ci(spoof_ics)
    print lb, given_ic, ub
    assert_less_equal(lb - 1, given_ic)
    assert_less_equal(given_ic, ub + 1)

def test_spoof_motif_cftp_occ():
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
    bio_matrix = matrix_from_motif(motif)
    N = len(motif)
    print occupancies
    mean_bio_occ = mean(occupancies(motif))
    motifs = spoof_motif_cftp_occ(motif, num_motifs=100)
    spoof_occs = map(lambda m:mean(occupancies(m)), motifs)
    lb, ub = mean_ci(spoof_occs)
    print lb, mean_bio_occ, ub
    assert_less_equal(lb, mean_bio_occ)
    assert_less_equal(mean_bio_occ, ub)

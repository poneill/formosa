from maxent_sampling import spoof_maxent_motifs
from formosa_utils import motif_mi, count, mi_test_cols, transpose, choose2
import cPickle
import random
from utils import fdr
import sys
from tqdm import *

def mi_percentile(motif, trials=1000, verbose=False):
    obs_mi = motif_mi(motif)
    iterator = tqdm if verbose else lambda x:x
    spoof_mis = map(motif_mi, iterator(spoof_maxent_motifs(motif,num_motifs=trials)))
    return count(lambda x:x >= obs_mi, spoof_mis)/float(trials)

def motif_mi_col_test(motif, trials=1000):
    cols = transpose(motif)
    return sum(mi_test_cols(colA, colB) for colA, colB in choose2(cols))/float(len(choose2(cols)))

def motif_test_cols(motif):
    cols = transpose(motif)
    return [mi_test_cols(colA, colB, alpha=None) for colA, colB in choose2(cols)]
    
def motif_mi_distances(motif, trials=1000):
    cols = transpose(motif)
    L = len(cols)
    correlated_distances = [j-i for (i,coli), (j,colj) in choose2(list(enumerate(cols)))
                            if mi_test_cols(coli, colj)]
    return (correlated_distances, L)
    
def get_motifs():
    sys.path.append("/home/pat/jaspar")
    from parse_jaspar import jaspar_motifs as euk_motifs
    sys.path.append("/home/pat/motifs")
    from parse_tfbs_data import tfdf
    prok_motifs = [getattr(tfdf,tf) for tf in tfdf.tfs]
    return prok_motifs, euk_motifs

def do_mi_tests():
    random.seed("do_mi_tests")
    prok_motifs, euk_motifs = get_motifs()
    prok_tests = [motif_test_cols(motif) for motif in tqdm(prok_motifs)]
    euk_tests = [motif_test_cols(motif) for motif in tqdm(euk_motifs)]
    with open("prok_tests.pkl",'w') as f:
        cPickle.dump(prok_tests,f)
    with open("euk_tests.pkl",'w') as f:
        cPickle.dump(euk_tests,f)

def analyze_mi_tests(prok_tests, euk_tests):
    pass
    prok_q = fdr(concat(prok_tests))
    euk_q = fdr(concat(euk_tests))
    prok_correlated_percentage = count(lambda x:x <= q,(concat(prok_tests)))/float(len(concat(prok_tests)))
    euk_correlated_percentage = count(lambda x:x <= q,(concat(euk_tests)))/float(len(concat(euk_tests)))
    prok_ds = [[j - i for (i, coli), (j,colj) in choose2(list(enumerate(transpose(motif))))]
               for motif in prok_motifs]
    euk_ds = [[j - i for (i, coli), (j,colj) in choose2(list(enumerate(transpose(motif))))]
               for motif in euk_motifs]
    def binom_ci(xs):
        """return width of error bar"""
        bs_means = sorted([mean(bs(xs)) for x in range(1000)])
        mu = mean(xs)
        return (mu - bs_means[25], bs_means[975] - mu)
    prok_cis = [binom_ci([t <= prok_q for t,d in zip(concat(prok_tests), concat(prok_ds)) if d == i])
                for i in trange(1,20)]
    euk_cis = [binom_ci([t <= euk_q for t,d in zip(concat(euk_tests), concat(euk_ds)) if d == i])
                for i in trange(1,20)]
    plt.errorbar(range(1,20),[mean([t <= prok_q for t,d in zip(concat(prok_tests), concat(prok_ds)) if d == i]) for i in range(1,20)],yerr=transpose(prok_cis),label="Prokaryotic Motifs")
    plt.errorbar(range(1,20),[mean([t <= euk_q for t,d in zip(concat(euk_tests), concat(euk_ds)) if d == i]) for i in range(1,20)],yerr=transpose(euk_cis),label="Eukaryotic Motifs")
    plt.xlabel("Distance (bp)",fontsize="large")
    plt.ylabel("Proportion of Significant Correlations",fontsize="large")
    plt.legend()
    
    
    

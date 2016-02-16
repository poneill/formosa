from maxent_sampling import spoof_maxent_motifs
from formosa_utils import motif_mi, count, mi_test_cols, transpose, choose2
import cPickle
import random
from utils import fdr, inverse_cdf_sample, normalize, concat, maybesave
import sys
from tqdm import *
from math import sqrt
from matplotlib import pyplot as plt
from collections import defaultdict
import seaborn as sns

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

def get_tests():
    with open("prok_tests.pkl") as f:
        prok_tests = cPickle.load(f)
    with open("euk_tests.pkl") as f:
        euk_tests = cPickle.load(f)
    return prok_tests, euk_tests
        
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
    plt.errorbar(range(1,20),
                 [mean([t <= prok_q for t,d in zip(concat(prok_tests), concat(prok_ds)) if d == i])
                  for i in range(1,20)],yerr=transpose(prok_cis),label="Prokaryotic Motifs",capthick=1)
    plt.errorbar(range(1,20),
                 [mean([t <= euk_q for t,d in zip(concat(euk_tests), concat(euk_ds)) if d == i])
                  for i in range(1,20)],yerr=transpose(euk_cis),label="Eukaryotic Motifs",capthick=1)
    plt.xlabel("Distance (bp)",fontsize="large")
    plt.ylabel("Proportion of Significant Correlations",fontsize="large")
    plt.legend()

def analyze_correlated_digrams(prok_tests, euk_tests, filename=None):
    digrams = [(b1,b2) for b1 in "ACGT" for b2 in "ACGT"]
    prok_q = fdr(concat(prok_tests))
    euk_q = fdr(concat(euk_tests))
    prok_digrams = defaultdict(int)
    prok_corr_digrams = defaultdict(int)
    prok_adj_digrams = defaultdict(int)
    for tests, motif in tqdm(zip(prok_tests, prok_motifs)):
        for test, ((i,coli),(j,colj)) in zip(tests, choose2(list(enumerate(transpose((motif)))))):
            for bi,bj in transpose((coli,colj)):
                prok_digrams[(bi,bj)] += 1
                if j == i + 1:
                    prok_adj_digrams[(bi,bj)] += 1
                if test <= prok_q:
                    prok_corr_digrams[(bi,bj)] += 1
    prok_corr_N = float(sum(prok_corr_digrams.values()))
    prok_adj_N = float(sum(prok_adj_digrams.values()))
    prok_N = float(sum(prok_digrams.values()))
    prok_ps = normalize(prok_digrams.values())
    prok_adj_ps = normalize(prok_adj_digrams.values())
    prok_corr_ps = normalize(prok_corr_digrams.values())
    prok_yerr = [1.96*sqrt(1.0/prok_N*p*(1-p)) for p in prok_ps]
    prok_adj_yerr = [1.96*sqrt(1.0/prok_adj_N*p*(1-p)) for p in prok_adj_ps]
    prok_corr_yerr = [1.96*sqrt(1.0/prok_corr_N*p*(1-p)) for p in prok_corr_ps]

    euk_digrams = defaultdict(int)
    euk_corr_digrams = defaultdict(int)
    euk_adj_digrams = defaultdict(int)
    for tests, motif in tqdm(zip(euk_tests, euk_motifs)):
        for test, ((i,coli),(j,colj)) in zip(tests, choose2(list(enumerate(transpose((motif)))))):
            for bi,bj in transpose((coli,colj)):
                euk_digrams[(bi,bj)] += 1
                if j == i + 1:
                    euk_adj_digrams[(bi,bj)] += 1
                if test <= euk_q:
                    euk_corr_digrams[(bi,bj)] += 1
    euk_corr_N = float(sum(euk_corr_digrams.values()))
    euk_adj_N = float(sum(euk_adj_digrams.values()))
    euk_N = float(sum(euk_digrams.values()))
    euk_ps = normalize(euk_digrams.values())
    euk_adj_ps = normalize(euk_adj_digrams.values())
    euk_corr_ps = normalize(euk_corr_digrams.values())
    euk_yerr = [1.96*sqrt(1.0/euk_N*p*(1-p)) for p in euk_ps]
    euk_adj_yerr = [1.96*sqrt(1.0/euk_adj_N*p*(1-p)) for p in euk_adj_ps]
    euk_corr_yerr = [1.96*sqrt(1.0/euk_corr_N*p*(1-p)) for p in euk_corr_ps]

    palette = sns.cubehelix_palette(4)
    ax = plt.subplot(211)
    # plt.bar(range(16),normalize(prok_digrams.values()))
    # plt.bar(range(16),normalize(prok_corr_digrams.values()),color='g')
    # plt.bar([x-0.2 for x in range(16)], prok_relative_ratios.values(), color='g', label="Correlated Column-pairs",width=0.2)
    # plt.bar([x for x in range(16)],prok_adj_relative_ratios.values(),color='r',alpha=1,yerr=prok_adj_yerr,label="Adjacent Column-pairs",width=0.2)
    # plt.bar([x+0.2 for x in range(16)],[1]*16,color='b',alpha=1,yerr=(prok_yerr),capsize=10,capstyle='butt',label="All Column-pairs",width=0.2)
    plt.bar([x-0.2 for x in range(16)], prok_ps, label="All Column-Pairs",width=0.2,yerr=prok_yerr,color=palette[0])
    plt.bar([x for x in range(16)],normalize(prok_adj_digrams.values()),label="Adjacent Column-Pairs",
            width=0.2,yerr=prok_adj_yerr,color=palette[1])
    plt.bar([x+0.2 for x in range(16)],normalize(prok_corr_digrams.values()),alpha=1,
            capstyle='butt',label="Correlated Column-Pairs",width=0.2,yerr=prok_corr_yerr,color=palette[3])
    plt.plot([0,16],[1.0/16, 1.0/16],linestyle='--',color=palette[3],label="Equiprobability",linewidth=1)
    ax.set_xticks([x for x in range(16)])
    ax.set_xticklabels( ["".join(dg) for dg in digrams],fontsize='large')
    plt.xlim(-0.5,15.5)
    plt.ylim(0,0.2)
    #plt.xlabel("Dimer",fontsize='large')
    plt.ylabel("Prokaryotic Frequency",fontsize='large')
    #plt.ylim(0,2)
    plt.legend()
    
    ax2 = plt.subplot(212)
    plt.plot([0,16],[1.0/16, 1.0/16],linestyle='--',color=palette[3],label="Equiprobability",linewidth=1)
    plt.bar([x-0.2 for x in range(16)], euk_ps, label="All Column-Pairs",width=0.2,yerr=euk_yerr,color=palette[0])
    plt.bar([x for x in range(16)],normalize(euk_adj_digrams.values()),label="Adjacent Column-Pairs",
            width=0.2,yerr=euk_adj_yerr,color=palette[1])
    plt.bar([x+0.2 for x in range(16)],normalize(euk_corr_digrams.values()),alpha=1,
            capstyle='butt',label="Correlated Column-Pairs",width=0.2,yerr=euk_corr_yerr,color=palette[3])
    ax2.set_xticks([x for x in range(16)])
    ax2.set_xticklabels( ["".join(dg) for dg in digrams],fontsize='large')
    #plt.xlabel("Dimer",fontsize='large')
    plt.xlim(-0.5,15.5)
    plt.ylim(0,0.2)
    plt.ylabel("Eukaryotic Frequency",fontsize='large')
    #plt.ylim(0,2)
    plt.legend()
    maybesave(filename)
    

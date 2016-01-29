from formosa_utils import h, log_fac, np_log_normalize, secant_interval
from formosa_utils import motif_ic, mean_ci, inverse_cdf_sampler, concat
from formosa_utils import permute, transpose
import numpy as np
from math import log, exp
from collections import Counter
from tqdm import *
import random

def maxent_motif(N,L,desired_ic,tolerance=10**-10,beta=None,verbose=False):
    """sample motif from max ent distribution with mean desired_ic"""
    # first we adjust the desired ic upwards so that when motif_ic is
    # called with 1st order correction, we get the desired ic.
    if beta is None:
        if verbose:
            print "finding beta"
        correction_per_col = 3/(2*log(2)*N)
        desired_ic += L * correction_per_col
        beta = find_beta_for_mean_motif_ic(N,L,desired_ic,tolerance=tolerance,verbose=verbose)
    ps = count_ps_from_beta(N,beta)
    count_sampler = inverse_cdf_sampler(enumerate_counts(N),ps)
    counts = [count_sampler() for i in range(L)]
    cols = [sample_col_from_count(count) for count in counts]
    return map(lambda site:"".join(site),transpose(cols))

def spoof_maxent_motif(motif,verbose=False):
    N = len(motif)
    L = len(motif[0])
    des_ic = motif_ic(motif)
    if verbose:
        print "n: {} L: {} des_ic: {}".format(n,L,des_ic)
    return maxent_motif(N,L,des_ic,verbose=verbose)

def spoof_maxent_motifs(motif,num_motifs,verbose=False):
    n = len(motif)
    L = len(motif[0])
    des_ic = motif_ic(motif)
    if verbose:
        print "n: {} L: {} des_ic: {}".format(n,L,des_ic)
    return maxent_motifs(n,L,des_ic,num_motifs,verbose=verbose)

def maxent_motifs(N, L, desired_ic, num_motifs, tolerance=10**-10,beta=None,verbose=False):
    if beta is None:
        correction_per_col = 3/(2*log(2)*N)
        desired_ic += L * correction_per_col
        beta = find_beta_for_mean_motif_ic(N,L,desired_ic,tolerance=tolerance,verbose=verbose)
        if verbose:
            print "beta:",beta
    ps = count_ps_from_beta(N,beta)
    count_sampler = inverse_cdf_sampler(enumerate_counts(N),ps)
    def sample():
        counts = [count_sampler() for i in range(L)]
        cols = [sample_col_from_count(count) for count in counts]
        return map(lambda site:"".join(site),transpose(cols))
    return [sample() for _ in trange(num_motifs)]
    
def find_beta_for_mean_motif_ic(n,L,desired_ic,tolerance=10**-10,verbose=False):
    desired_ic_per_col = desired_ic/L
    return find_beta_for_mean_col_ic(n,desired_ic_per_col,tolerance,verbose=verbose)

def find_beta_for_mean_col_ic(n, desired_ic_per_col,tolerance=10**-10,verbose=False):
    """find beta such that entropy*exp(-beta*entropy)/Z = des_ent"""
    if verbose:
        print "enumerating countses"
    countses = enumerate_counts(n)
    if verbose:
        print "enumerating entropies"
    entropies = np.array(map(entropy_from_counts, countses))
    #cols = np.array(map(countses_to_cols, countses))
    if verbose:
        print "enumerating cols"
    #cols = np.exp(np.array(map(log_counts_to_cols, countses)))
    iterator = tqdm(countses) if verbose else countses
    log_cols = np.array(map(log_counts_to_cols, iterator))
    def f2(beta):
        log_phats = np_log_normalize(log_cols + -beta*entropies)
        expected_entropy = np.exp(log_phats).dot(entropies)
        return 2 - expected_entropy - desired_ic_per_col
    ub = 1000
    while f2(ub) < 0:
        ub *= 2
        print "raising upper bound to:",ub
    return secant_interval(f2,0,ub,verbose=verbose,tolerance=tolerance)

def enumerate_counts(N):
    return list(partitionfunc(N,4,l=0))

def enumerate_counts_iter(N):
    return (partitionfunc(N,4,l=0))
    
def partitionfunc(n,k,l=1):
    '''n is the integer to partition, k is the length of partitions, l is the min partition element size'''
    # borrowed from stack overflow
    if k < 1:
        raise StopIteration
    if k == 1:
        if n >= l:
            yield (n,)
        raise StopIteration
    for i in range(l,n+1):
        for result in partitionfunc(n-i,k-1,i):
            yield (i,)+result

def count_ps_from_beta(n,beta,verbose=True):
    log_ws = np.array([log_counts_to_cols(count) + (-beta*entropy_from_counts(count))
              for count in tqdm(enumerate_counts_iter(n))])
    return np.exp(np_log_normalize(log_ws))

def entropy_from_counts(counts):
    N = float(sum(counts))
    ps = [c/N for c in counts]
    return h(ps)

def log_counts_to_cols(counts):
    """return number of cols associated given counts"""
    N = sum(counts)
    # all_cols = 4**N
    metacounts = Counter(counts)
    log_counts_to_bases = log_fac(4) - sum(log_fac(multiplicity) for multiplicity in metacounts.values())
    log_bases_to_pos = (log_fac(N) - sum(log_fac(count) for count in counts))
    return log_counts_to_bases + log_bases_to_pos

def sample_col_from_count(count):
    col = concat([[base]*n for base,n in zip(permute("ACGT"),count)])
    random.shuffle(col)
    return col

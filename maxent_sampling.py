from formosa_utils import h, log_fac, np_log_normalize, secant_interval
from formosa_utils import motif_ic, mean_ci, inverse_cdf_sampler, concat
from formosa_utils import permute, transpose
import numpy as np
from math import log, exp
from collections import Counter, defaultdict
from tqdm import *
import random

log2 = lambda x:log(x,2)
dna_alphabet = "ACGT"
aa_alphabet = "ACDEFGHIKLMNPQRSTVWY"

def logger(s, verbose):
    if verbose:
        print s
        
def maxent_motif(N,L,desired_ic,tolerance=10**-10,beta=None,verbose=False, A=4):
    """sample motif from max ent distribution with mean desired_ic"""
    # first we adjust the desired ic upwards so that when motif_ic is
    # called with 1st order correction, we get the desired ic.
    if beta is None:
        if verbose:
            print "finding beta"
        correction_per_col = (A-1)/(2*log(2)*N)
        desired_ic += L * correction_per_col
        beta = find_beta_for_mean_motif_ic(N,L,desired_ic,tolerance=tolerance,verbose=verbose, A=A)
    ps = count_ps_from_beta(N, beta, A=A)
    count_sampler = inverse_cdf_sampler(enumerate_counts(N, A), ps)
    counts = [count_sampler() for i in range(L)]
    cols = [sample_col_from_count(count, A=A) for count in counts]
    return map(lambda site:"".join(site),transpose(cols))

def is_dna(xs):
    if type(xs) is str:
        return all(c in dna_alphabet for c in xs)
    elif type(xs) is list:
        return all(is_dna(x) for x in xs)
        
def spoof_maxent_motif(motif,verbose=False, A=4):
    N = len(motif)
    L = len(motif[0])
    des_ic = motif_ic(motif)
    if verbose:
        print "n: {} L: {} des_ic: {}".format(n,L,des_ic)
    return maxent_motif(N,L,des_ic,verbose=verbose, A=A)

def spoof_maxent_motifs(motif,num_motifs,verbose=False, A=4):
    N = len(motif)
    L = len(motif[0])
    des_ic = motif_ic(motif, A=A)
    if verbose:
        print "N: {} L: {} des_ic: {}".format(N,L,des_ic)
    return maxent_motifs(N,L,des_ic,num_motifs,verbose=verbose, A=A)

def maxent_motifs(N, L, desired_ic, num_motifs, tolerance=10**-10, A=4, beta=None, countses=None,
                  entropies=None, log_cols=None, verbose=False):
    ### computational
    if countses is None:
        logger("countses", verbose)
        countses = enumerate_counts(N, A, verbose=verbose)
    if entropies is None:
        logger("entropies", verbose)
        entropies = np.array(map(entropy_from_counts, tqdm(countses)))
    if log_cols is None:
        iterator = tqdm(countses) if verbose else countses
        logger("log_cols", verbose)
        log_cols = np.array([log_counts_to_cols(counts, A=A) for counts in iterator])
    if beta is None:
        correction_per_col = (A-1)/(2*log(2)*N)
        desired_ic += L * correction_per_col
        beta = find_beta_for_mean_motif_ic(N,L,desired_ic,tolerance=tolerance,verbose=verbose, A=A,
                                           countses=countses, entropies=entropies, log_cols=log_cols)
        logger("beta: %s" % beta, verbose)
    logger("computing count ps from beta", verbose)
    ps = count_ps_from_beta(N,beta, A=A, verbose=verbose,
                            log_cols=log_cols, entropies=entropies)
    count_sampler = inverse_cdf_sampler(countses, ps)
    def sample():
        counts = [count_sampler() for i in range(L)]
        cols = [sample_col_from_count(count, A=A) for count in counts]
        return map(lambda site:"".join(site),transpose(cols))
    iterator = trange if verbose else xrange
    if verbose:
        print "sampling"
    return [sample() for _ in iterator(num_motifs)]
    
def find_beta_for_mean_motif_ic(n,L,desired_ic,tolerance=10**-10,verbose=False, A=4,
                                countses=None, entropies=None, log_cols=None):
    desired_ic_per_col = desired_ic/float(L)
    return find_beta_for_mean_col_ic(n,desired_ic_per_col,tolerance,verbose=verbose, A=A,
                                     countses=countses, entropies=entropies, log_cols=log_cols)

def find_beta_for_mean_col_ic(n, desired_ic_per_col,tolerance=10**-10,verbose=False, A=4,
                              countses=None, entropies=None, log_cols=None):
    """find beta such that entropy*exp(-beta*entropy)/Z = des_ent"""
    if countses is None:
        if verbose:
            print "enumerating countses"
        countses = enumerate_counts(n, A, verbose=verbose)
    if entropies is None:
        if verbose:
            print "enumerating entropies"
            entropies = np.array(map(entropy_from_counts, tqdm(countses)))
        else:
            entropies = np.array(map(entropy_from_counts, countses))
    #cols = np.array(map(countses_to_cols, countses))
    if log_cols is None:
        if verbose:
            print "enumerating cols"
            #cols = np.exp(np.array(map(log_counts_to_cols, countses)))
        iterator = tqdm(countses) if verbose else countses
        log_cols = np.array([log_counts_to_cols(counts, A=A) for counts in iterator])
    def f2(beta):
        log_phats = np_log_normalize(log_cols + -beta*entropies)
        expected_entropy = np.exp(log_phats).dot(entropies)
        return log2(A) - expected_entropy - desired_ic_per_col
    lb = -1
    while f2(lb) > 0:
        lb *= 2
        if lb < -1000:
            print "Warning, failed to find lower bound on beta"
            raise Exception("Couldn't find beta'")
    ub = 1000
    while f2(ub) < 0:
        ub *= 2
        print "raising upper bound to:",ub
    return secant_interval(f2,lb,ub,verbose=verbose,tolerance=tolerance)

def enumerate_counts(N, A=4, verbose=False):
    if verbose:
        return list(tqdm(partitionfunc(N,A,l=0), total=partitions(N, A)))
    else:
        return list(partitionfunc(N,A,l=0))

def exact_partitions(N, k):
    """compute number of partitions of N into exactly k parts"""
    parts = {(0,0):1}
    def p(N, k):
        #print "starting on %s, %s" % (N, k)
        if (N, k) in parts:
            #print "looking up"
            ans = parts[N, k]
        elif N <= 0 or k <= 0:
            ans = 0
        else:
            ans = p(N-k,k) + p(N-1, k-1)
            parts[N, k] = ans
        #print "%s, %s: %s" % (N, k, ans)
        return ans
    return p(N, k)
            
def partitions(N, k):
    """compute number of partitions of N into at most k parts"""
    return sum(exact_partitions(N, j) for j in range(1, k+1))
    
def enumerate_aa_counts(N):
    return list(partitionfunc(N,20,l=0))
    
def enumerate_counts_iter(N, A=4, verbose=False):
    if verbose:
        return tqdm(partitionfunc(N,A,l=0), total=partitions(N, A))
    else:
        return (partitionfunc(N,A,l=0))
    
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

def count_ps_from_beta(N, beta, verbose=False, A=4, log_cols=None, entropies=None):
    iterator = (lambda x:tqdm(x, total=partitions(N, A))) if verbose else (lambda x:x)
    if log_cols is None or entropies is None:
        log_ws = np.array([log_counts_to_cols(count, A=A) + (-beta*entropy_from_counts(count))
                           for count in iterator(enumerate_counts_iter(N, A=A))])
    else:
        log_ws = log_cols + -beta*entropies
    return np.exp(np_log_normalize(log_ws))

def entropy_from_counts(counts):
    N = float(sum(counts))
    ps = [c/N for c in counts]
    return h(ps)

def log_counts_to_cols_ref(counts, A=4):
    """return number of cols associated given counts"""
    N = sum(counts)
    # all_cols = 4**N
    metacounts = Counter(counts)
    log_counts_to_bases = log_fac(A) - sum(log_fac(multiplicity) for multiplicity in metacounts.values())
    log_bases_to_pos = (log_fac(N) - sum(log_fac(count) for count in counts))
    return log_counts_to_bases + log_bases_to_pos

def log_counts_to_cols(counts, A=4):
    N = sum(counts)
    # all_cols = 4**N
    metacounts = defaultdict(int)
    for c in counts:
        metacounts[c] += 1
    log_counts_to_bases = log_fac(A) - sum(log_fac(multiplicity) for multiplicity in metacounts.values())
    log_bases_to_pos = (log_fac(N) - sum(log_fac(count) for count in counts))
    return log_counts_to_bases + log_bases_to_pos


def test_log_counts_to_cols(A=4):
    counts = sorted(tuple(random.randrange(5) for i in range(A)))
    assert log_counts_to_cols(counts, A) == log_counts_to_cols_ref(counts, A)

def sample_col_from_count(count, A=4):
    if A == 4:
        alphabet = dna_alphabet
    elif A == 20:
        alphabet = aa_alphabet
    col = concat([[base]*n for base,n in zip(permute(alphabet),count)])
    random.shuffle(col)
    return col

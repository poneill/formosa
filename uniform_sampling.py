from formosa_utils import transpose, inverse_cdf_sampler, sample_until, inrange
from formosa_utils import motif_ic
from maxent_sampling import sample_col_from_count, find_beta_for_mean_motif_ic
from maxent_sampling import count_ps_from_beta, enumerate_counts
from math import exp, log
from tqdm import *
import random

def uniform_motif(N,L,desired_ic,epsilon=0.1,beta=None,ps=None,count_sampler=None,verbose=False):
    if verbose:  print "uniform motif accept reject:",N,L,desired_ic,beta
    correction_per_col = 3/(2*log(2)*N)
    desired_ic_for_beta = desired_ic + L * correction_per_col
    if desired_ic_for_beta == 2*L: # if we reach the upper limit, things break down
        cols = [sample_col_from_count((0,0,0,N)) for _ in range(L)]
        motif_p = map(lambda site:"".join(site),transpose(cols))
        return motif_p
    if beta is None:
        beta = find_beta_for_mean_motif_ic(N,L,desired_ic_for_beta)
        if verbose:
            print "beta:",beta
    if ps is None:
        ps = count_ps_from_beta(N,beta)
    if count_sampler is None:
        count_sampler = inverse_cdf_sampler(enumerate_counts(N),ps)
    def rQ_raw():
        counts = [count_sampler() for i in range(L)]
        cols = [sample_col_from_count(count) for count in counts]
        motif_p = map(lambda site:"".join(site),transpose(cols))
        return motif_p
    def rQ():
        return sample_until(lambda M:inrange(M,desired_ic,epsilon),rQ_raw,1,progress_bar=False)[0]
    def dQhat(motif):
        return exp(beta*motif_ic(motif))
    Imin = desired_ic - epsilon
    Imax = desired_ic + epsilon
    log_M = -beta*Imin
    if verbose: print "Imin, Imax, log_M:",Imin, Imax, log_M
    def dQ(motif):
        return exp(beta*motif_ic(motif) + log_M)
    def AR(motif):
        return 1.0/dQ(motif)
    #M = exp(-beta*(desired_ic - epsilon)) # which ic? +/- correction
    trials = 0
    while True:
        trials +=1
        motif = rQ()
        r = random.random()
        if r < AR(motif):
            return motif
        if verbose and trials % 100 == 0:
            print trials, AR(motif)

def uniform_motifs(N,L,desired_ic,num_motifs,epsilon=0.1,beta=None,verbose=False):
    if beta is None:
        correction_per_col = 3/(2*log(2)*N)
        desired_ic_for_beta = desired_ic + L * correction_per_col
        beta = find_beta_for_mean_motif_ic(N,L,desired_ic_for_beta,verbose=verbose)
    ps = count_ps_from_beta(N,beta)
    count_sampler = inverse_cdf_sampler(enumerate_counts(N),ps)
    return [uniform_motif(N,L,desired_ic,epsilon=epsilon,beta=beta,
                                        ps=ps,count_sampler=count_sampler,verbose=verbose)
            for i in trange(num_motifs)]

def spoof_uniform_motifs(motif, num_motifs, epsilon=0.1,verbose=False):
    N, L = len(motif), len(motif[0])
    desired_ic = motif_ic(motif)
    if verbose: print "starting spoof motifs uniform with:",N, L, desired_ic
    return uniform_motifs(N, L, desired_ic, num_motifs, epsilon, verbose=verbose)

def spoof_uniform_motif(motif, epsilon=0.1, verbose=False):
    N, L = len(motif), len(motif[0])
    desired_ic = motif_ic(motif)
    if verbose: print "starting spoof motifs uniform with:",n,L,desired_ic
    return uniform_motifs(N, L, desired_ic, num_motifs, epsilon, verbose=verbose)[0]

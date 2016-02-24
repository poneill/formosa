"""perform exact evolutionary sampling via CFTP"""

import random
from tqdm import *
from math import exp, log
from scipy import polyfit, poly1d
from formosa_utils import sorted_indices, subst, mutate_site, seq_scorer, argmin
from formosa_utils import argmax, rslice, sigma_from_matrix, pssm_from_motif
from formosa_utils import psfm_from_motif, motif_ic, sample_matrix, approx_mu
from formosa_utils import mean, matrix_from_motif, occupancies

def sample_motif_cftp(matrix, mu, Ne, n,verbose=False):
    iterator = trange(n,desc="sampling cftp motif") if verbose else xrange(n)
    return [sample_site_cftp(matrix, mu, Ne)
            for i in iterator]

def sample_site_cftp(matrix, mu, Ne):
    L = len(matrix)
    f = seq_scorer(matrix)
    def log_phat(s):
        ep = f(s)
        nu = Ne - 1
        return -nu*log(1 + exp(ep - mu))
    first_site = "A"*L
    last_site = "T"*L
    best_site = "".join(["ACGT"[argmin(row)] for row in matrix])
    worst_site = "".join(["ACGT"[argmax(row)] for row in matrix])
    #middle_sites  = [[random_site(L)] for i in range(10)]
    #trajs = [[best_site]] + middle_sites + [[worst_site]]
    trajs = [[best_site],[worst_site]]
    ords = [rslice("ACGT",sorted_indices(row)) for row in matrix]
    def mutate_site(site,(ri,direction)):
        b = (site[ri])
        idx = ords[ri].index(b)
        idxp = min(max(idx + direction,0),3)
        bp = ords[ri][idxp]
        return subst(site,bp,ri)
    iterations = 1
    rs = [(random.randrange(L),random.choice([-1,1]),random.random())
          for i in range(iterations)]
    converged = False
    while not converged:
        for ri, rdir, r in rs:
            for traj in trajs:
                x = traj[-1]
                xp = mutate_site(x,(ri, rdir))
                if log(r) < log_phat(xp) - log_phat(x):
                    x = xp
                traj.append(x)
        if trajs[0][-1] == trajs[-1][-1]:
            converged = True
        iterations *= 2
        #print iterations,[traj[-1] for traj in trajs]
        rs = [(random.randrange(L),random.choice([-1,1]),random.random())
              for i in range(iterations)] + rs
    assert all(map(lambda traj:traj[-1] == trajs[0][-1],trajs))
    return trajs[0][-1]
    #return trajs

def spoof_motif_cftp(motif, num_motifs=10, trials=1, sigma=None,Ne_tol=10**-2,verbose=False):
    n = len(motif)
    L = len(motif[0])
    copies = 10*n
    if sigma is None: sigma = sigma_from_matrix(pssm_from_motif(motif,pc=1))
    print "sigma:", sigma
    bio_ic = motif_ic(motif)
    matrix = sample_matrix(L, sigma)
    mu = approx_mu(matrix, copies=10*n, G=5*10**6)
    print "mu:", mu
    def f(Ne):
        motifs = [sample_motif_cftp(matrix, mu, Ne, n, verbose=verbose)
                  for i in trange(trials)]
        return mean(map(motif_ic,motifs)) - bio_ic
    # lb = 1
    # ub = 10
    # while f(ub) < 0:
    #     ub *= 2
    #     print ub
    x0s = [2,10]#(lb + ub)/2.0
    # print "choosing starting seed for Ne"
    # fs = map(lambda x:abs(f(x)),x0s)
    # print "starting values:",x0s,fs
    # x0 = x0s[argmin(fs)]
    # print "chose:",x0
    # Ne = bisect_interval_noisy_ref(f,x0,lb=1,verbose=True)
    Ne = log_regress_spec2(f,x0s,tol=Ne_tol)
    print "Ne:",Ne
    return [sample_motif_cftp(matrix, mu, Ne, n) for _ in trange(num_motifs)]

def spoof_motif_cftp_occ(motif, num_motifs=10, trials=1, sigma=None,Ne_tol=10**-2,verbose=False):
    """spoof motifs based on occupancy rather than motif IC"""
    N = len(motif)
    L = len(motif[0])
    copies = 10*N
    pssm = pssm_from_motif(motif,pc=1)
    if sigma is None: sigma = sigma_from_matrix(pssm)
    print "sigma:", sigma
    matrix = sample_matrix(L, sigma)
    bio_matrix = matrix_from_motif(motif)
    mu = approx_mu(matrix, copies=copies, G=5*10**6)
    mean_bio_occ = mean(occupancies(motif))
    print "mu:", mu
    def f(Ne):
        motifs = [sample_motif_cftp(matrix, mu, Ne, N, verbose=verbose)
                  for i in trange(trials)]
        
        return mean(map(lambda m:mean(occupancies(m)), motifs)) - mean_bio_occ
    # lb = 1
    # ub = 10
    # while f(ub) < 0:
    #     ub *= 2
    #     print ub
    x0s = [2,10]#(lb + ub)/2.0
    # print "choosing starting seed for Ne"
    # fs = map(lambda x:abs(f(x)),x0s)
    # print "starting values:",x0s,fs
    # x0 = x0s[argmin(fs)]
    # print "chose:",x0
    # Ne = bisect_interval_noisy_ref(f,x0,lb=1,verbose=True)
    Ne = log_regress_spec2(f,x0s,tol=Ne_tol)
    print "Ne:",Ne
    return [sample_motif_cftp(matrix, mu, Ne, N) for _ in trange(num_motifs)]

def log_regress_spec2(f,xs, tol=0.01, diagnose=False):
    """find root f(x) = 0 using logistic regression, starting with xs, using weighted regression"""
    print "initial seeding for log_regress (weighted)"
    ys = map(f,xs)
    log_xs = map(log,xs)
    plotting = False
    honest_guesses = []
    while len(honest_guesses) < 2 or abs(honest_guesses[-1] -
                                     honest_guesses[-2]) > tol:
        m, b = weighted_regress(log_xs,ys)
        honest_guess = -b/m
        dx = 0#-(honest_guesses[-1] - honest_guess) if honest_guesses else 0
        log_xp = honest_guess + 2*dx
        log_xs.append(log_xp)
        yxp = f(exp(log_xp))
        ys.append(yxp)
        honest_guesses.append(honest_guess)
        diff = (abs((honest_guesses[-1]) - (honest_guesses[-2]))
                if len(honest_guesses) > 1 else None)
        print "honest_guess:",(honest_guess),"xp:",(log_xp),\
            "y:",yxp, "diff:",diff
    m, b = weighted_regress(log_xs,ys)
    log_xp = -b/m
    print "final guess: log_xp:",log_xp
    if diagnose:
        return log_xs,ys
    else:
        return exp(log_xp)

def weighted_regress(xs,ys):
    avg_yval = mean(map(abs,ys))
    ws = [exp(-abs(y)/avg_yval) for y in ys]
    return (polyfit(xs,ys,1,w=ws))

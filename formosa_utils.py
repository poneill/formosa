from math import log, sqrt
from scipy.special import gammaln
import numpy as np
import bisect
import random

def h(ps):
    """compute entropy (in bits) of a probability distribution ps"""
    return -sum([p * safe_log2(p) for p in ps])

def safe_log2(x):
    """Implements log2, but defines log2(0) = 0"""
    return log(x,2) if x > 0 else 0

def log_fac(n):
    return gammaln(n+1)

def np_log_normalize(log_xs):
    log_Z = np_log_sum(log_xs)
    return log_xs - log_Z

def np_log_sum(log_xs):
    "given numpy array log_xs, return log(sum(xs))"
    log_xmax = np.max(log_xs)
    return log_xmax + log(np.sum(np.exp(log_xs - log_xmax)))

def secant_interval(f,xmin,xmax,ymin=None,ymax=None,tolerance=1e-10,verbose=False):
    if ymin is None:
        ymin = f(xmin)
    if ymax is None:
        ymax = f(xmax)
    #for iteration in xrange(1000):
    y = 1000000
    while abs(y) > tolerance:
        if verbose:
            print xmin,xmax,ymin,ymax
        assert(np.sign(ymin)!= np.sign(ymax)), "ymin=%s,ymax=%s" % (ymin,ymax)
        m = (ymax - ymin)/float(xmax - xmin)
        x = xmax - ymax/m
        y = f(x)
        if abs(y) < tolerance:
            return x
        else:
            if np.sign(y) == np.sign(ymin):
                xmin = x
                ymin = y
            else:
                xmax = x
                ymax = y

def inverse_cdf_sampler(xs,ps):
    """make a bintree for Sampling from discrete distribution ps over set xs"""
    cum_ps = cumsum(ps)
    def sampler():
        r = random.random()
        i = bisect.bisect_left(cum_ps,r)
        return xs[i]
    return sampler

def cumsum(xs):
    acc = 0
    acc_list = []
    for x in xs:
        acc += x
        acc_list.append(acc)
    return acc_list

def concat(xss):
    return sum(xss,[])
    #return [x for xs in xxs for x in xs]

def permute(xs):
    """Return a random permutation of xs"""
    #return np.random.permutation(xs)
    xs_ = list(xs[:])
    random.shuffle(xs_)
    return xs_

def transpose(xxs):
    return zip(*xxs)

def motif_ic(motif,correct=True,alphabet_size=4):
    """Return the entropy of a motif, assuming independence and a
    uniform genomic background"""
    site_length = len(motif[0])
    return (log2(alphabet_size) * site_length -
            motif_entropy(motif,correct=correct,alphabet_size=4))

def motif_entropy(motif,correct=True,alphabet_size=4):
    """Return the entropy of a motif, assuming independence"""
    return sum(map(lambda col:entropy(col,correct=correct,alphabet_size=alphabet_size),
                   transpose(motif)))

def entropy(xs,correct=True,alphabet_size=None):
    """compute entropy (in bits) of a sample from a categorical
    probability distribution"""
    if alphabet_size == None:
        alphabet_size = len(set(xs)) # NB: assuming every element appears!
    ps = frequencies(xs)
    correction = ((alphabet_size - 1)/(2*log(2)*len(xs)) if correct
                  else 0) #Basharin 1959
    #print "correction:",correction
    return h(ps) + correction

def frequencies(xs):
    # faster than either of frequencies_ref [!]
    length = float(len(xs))
    return [xs.count(x)/length for x in set(xs)]
    
def log2(x):
    return log(x,2)

def mean_ci(xs):
    """Return 95% CI for mean"""
    mu = mean(xs)
    s = 1.96 * se(xs)
    return (mu - s,mu + s)

def mean(xs):
    if hasattr(xs,"__len__"):
        return sum(xs)/float(len(xs))
    else:
        acc = 0
        n = 0
        for x in xs:
            acc += x
            n += 1
        return acc/float(n)

def variance(xs,correct=True):
    n = len(xs)
    correction = n/float(n-1) if correct else 1
    mu = mean(xs)
    return correction * mean([(x-mu)**2 for x in xs])
    
def sd(xs,correct=True):
    return sqrt(variance(xs,correct=correct))
    
def se(xs,correct=True):
    return sd(xs,correct)/sqrt(len(xs))

def sample_until(p,sampler,n,progress_bar=True):
    """return n samples from sampler satisfying predicate p"""
    def gen():
        while True:
            x = sampler()
            if p(x):
                return x
    pb = trange if progress_bar else xrange #optional progress bar
    return [gen() for i in pb(n)]

def inrange(M,I,epsilon):
    return abs(motif_ic(M)-I) < epsilon

def motif_mi(motif):
    cols = transpose(motif)
    return sum([mi(col1,col2) for (col1,col2) in choose2(cols)])

def mi(xs,ys,correct=True):
    """Compute mutual information (in bits) of samples from two
    categorical probability distributions, using small sample size
    correction for each entropy computation.
    """
    hx  = entropy(xs,correct=correct,alphabet_size=4)
    hy  = entropy(ys,correct=correct,alphabet_size=4)
    hxy = entropy(zip(xs,ys),correct=correct,alphabet_size=16)
    return hx + hy - hxy

def mi_test_cols(xs, ys, alpha=0.05, trials=1000, perm_test=True):
    """using permutation test, do (1-sided) significance test for MI.  If
    permute, compare to random permutation, otherwise random site"""
    N = len(xs)
    obs_mi = mi(xs,ys)
    nullx = lambda : permute(xs) if perm_test else random_site(N)
    nully = lambda : permute(ys) if perm_test else random_site(N)
    perm_replicates = (mi(nullx(),nully()) for i in xrange(trials))
    perc = len(filter(lambda x:x >= obs_mi,perm_replicates))/float(trials)
    return perc < alpha

def count(p,xs):
    """Count number of xs satisfying p"""
    return len(filter(p,xs))

def choose2(xs, gen=False):
    """return list of choose(xs, 2) pairs, retaining ordering on xs"""
    if gen:
        return ((x1, x2) for i, x1 in enumerate(xs) for x2 in xs[i+1:])
    else:
        return [(x1, x2) for i, x1 in enumerate(xs) for x2 in xs[i+1:]]

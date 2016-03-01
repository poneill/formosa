from math import log, sqrt
from scipy.special import gammaln
import numpy as np
import bisect
import random
from math import pi, log, exp

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
    #return list(np.random.permutation(xs))
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

def motif_mi(motif, correct=False):
    cols = transpose(motif)
    return sum([mi(col1,col2, correct=correct) for (col1,col2) in choose2(cols)])

def mi(xs,ys,correct=True):
    """Compute mutual information (in bits) of samples from two
    categorical probability distributions, using small sample size
    correction for each entropy computation.
    """
    hx  = entropy(xs,correct=correct,alphabet_size=4)
    hy  = entropy(ys,correct=correct,alphabet_size=4)
    hxy = entropy(zip(xs,ys),correct=correct,alphabet_size=16)
    return hx + hy - hxy

def mi_test_cols(xs, ys, alpha=None, trials=1000):
    """using permutation test, do (1-sided) significance test for MI.  If
    permute, compare to random permutation, otherwise random site"""
    N = len(xs)
    obs_mi = mi(xs,ys)
    perm_replicates = (mi(permute(xs),permute(ys)) for i in xrange(trials))
    perc = len(filter(lambda x:x >= obs_mi,perm_replicates))/float(trials)
    if alpha:
        return perc < alpha
    else:
        return perc

def count(p,xs):
    """Count number of xs satisfying p"""
    return len(filter(p,xs))

def choose2(xs, gen=False):
    """return list of choose(xs, 2) pairs, retaining ordering on xs"""
    if gen:
        return ((x1, x2) for i, x1 in enumerate(xs) for x2 in xs[i+1:])
    else:
        return [(x1, x2) for i, x1 in enumerate(xs) for x2 in xs[i+1:]]

def sorted_indices(xs):
    """Return a list of indices that puts xs in sorted order.
    E.G.: sorted_indices([40,10,30,20]) => [1,3,2,0]"""
    return [i for (i,v) in sorted(enumerate(xs),key=lambda(i,v):v)]

def subst(xs,ys,i):
    """Substitute substring ys in xs, starting at i"""
    if not (type(ys) is list or type(ys) is str):
        ys = [ys]
    return xs[:i] + ys + xs[i+len(ys):]

def mutate_site(site):
    i = random.randrange(len(site))
    b = site[i]
    new_b = random.choice([c for c in "ACGT" if not c == b])
    return subst(site,new_b,i)

def score_seq(matrix,seq):
    #base_dict = {'A':0,'C':1,'G':2,'T':3}
    def base_dict(b):
        if b <= "C":
            if b == "A":
                return 0
            else:
                return 1
        elif b == "G":
            return 2
        else:
            return 3
    ans = 0
    for i in xrange(len(seq)):
        ans += matrix[i][base_dict(seq[i])]
    return ans

def seq_scorer(matrix):
    """accept matrix, return a function scoring sites"""
    # when score_seq JUST ISN'T FAST ENOUGH
    base_dicts = [{b:row[j] for j,b in enumerate("ACGT")} for row in matrix]
    def f(site):
        ans = 0
        for i in xrange(len(site)):
            ans += base_dicts[i][site[i]]
        return ans
    return f

def argmin(xs):
    i,x = min(enumerate(xs),key= lambda (i,x):x)
    return i

def argmax(xs):
    i,x = max(enumerate(xs),key= lambda (i,x):x)
    return i

def rslice(xs,js):
    return [xs[j] for j in js]

def sigma_from_matrix(matrix):
    """given a GLE matrix, estimate standard deviation of cell weights,
    correcting for bias of sd estimate.  See:
    https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
    """
    c = 2*sqrt(2/(3*pi))
    return mean(map(lambda x:sd(x,correct=True),matrix))/c

def pssm_from_motif(motif, pc=1):
    psfm = psfm_from_motif(motif, pc)
    return [[log(f/0.25,2) for f in row] for row in psfm]

def psfm_from_motif(motif,pc=1):
    n = float(len(motif))
    cols = transpose(motif)
    return [[(col.count(b) + pc)/(n+4*pc) for b in "ACGT"] for col in cols]

def normalize(xs):
    total = float(sum(xs))
    return [x/total for x in xs]

def sample_matrix(L,sigma):
    return [[random.gauss(0,sigma) for j in range(4)] for i in range(L)]

def approx_mu(matrix, copies, G=5*10**6):
    Zb = Zb_from_matrix(matrix, G)
    return log(copies) - log(Zb)

def Zb_from_matrix(matrix,G):
    """calculate partition function"""
    L = len(matrix)
    Zb_hat = prod(sum(exp(-ep) for ep in col) for col in matrix)/(4**L)
    return G * Zb_hat

def prod(xs):
    acc = 1
    for x in xs:
        acc *= x
    return acc

def occupancy(site, matrix, copies, G=5*10**6):
    mu = approx_mu(matrix, copies, G=G)
    ep = score_seq(matrix, site)
    return 1/(1+exp(ep-mu))

def occupancies(motif, copy_factor=10, G=5*10**6):
    matrix = matrix_from_motif(motif)
    copies = copy_factor * len(motif)
    mu = approx_mu(matrix, copies, G=G)
    eps = [score_seq(matrix, site) for site in motif]
    return [1/(1+exp(ep-mu)) for ep in eps]
        
def matrix_from_motif(seqs,pc=1):
    cols = transpose(seqs)
    N = float(len(seqs))
    raw_mat = [[-log((col.count(b)+pc)/(N+4*pc)) for b in "ACGT"] for col in cols]
    # now normalize each column by the average value
    avg = mean(map(mean,raw_mat))
    return [[x-avg for x in row] for row in raw_mat]

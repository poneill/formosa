"""Use maxent sampling to derive pvalues for motifs based on ic"""

from maxent_sampling import *
from utils import log_sum, random_motif, mean, random_site
from matplotlib import pyplot as plt
from math import log, exp, sqrt, erf, pi
from tqdm import *
from scipy.special import erfcx

def ic_log_pvalue(N, L, des_ic, verbose=False, trials=100, method="ub"):
    print des_ic
    correction_per_col = 3/(2*log(2)*N)
    K = L * correction_per_col # correction per motif
    ic_for_beta = des_ic + K
    tolerance = 10**-10
    beta = find_beta_for_mean_motif_ic(N,L,ic_for_beta,tolerance=tolerance,verbose=verbose) # correct val of beta
    countses = enumerate_counts(N)
    entropies = np.array(map(entropy_from_counts, countses))
    iterator = tqdm(countses) if verbose else countses
    log_cols = np.array(map(log_counts_to_cols, iterator))
    log_Zq = log_sum(log_cols + -beta*entropies)*L
    log_Zp = N*L*log(4)
    #log_prefactor = log_Zq - log_Zp + beta*2*L
    log_prefactor = log_Zq - log_Zp + beta*(2*L-K)
    if method == "UB":
        log_expectation_ub = (-beta*(des_ic))
        log_pval_ub = log_prefactor + log_expectation_ub
        return log_pval_ub - log(2)
    elif method == "analytic":
        mu, sigma = calc_params(N, L, beta)
        log_expectation = log(compute_expectation_spec(beta, mu, sigma))
        log_pval = log_prefactor + log_expectation
        return log_pval
    else:
        ms = maxent_motifs(N, L, des_ic, trials, beta=beta)
        ics = map(motif_ic, ms)
        print "des_ic, mean ics:", des_ic, mean(ics)
        log_expectation = log_sum([-beta*ic for ic in ics if ic > des_ic]) - log(trials) # Xxx loss of precision
        log_pval = log_prefactor + log_expectation
        return log_pval

def calc_params(N, L, beta):
    """compute ic sd of maxent motifs"""
    countses = enumerate_counts(N)
    entropies = np.array(map(entropy_from_counts, countses))
    iterator = tqdm(countses)
    log_cols = np.array(map(log_counts_to_cols, iterator))
    log_phats = np_log_normalize(log_cols + -beta*entropies)
    ps = np.exp(log_phats)
    correction_per_col = 3/(2*log(2)*N)
    K = L * correction_per_col # correction per motif
    mu = ps.dot(entropies)
    sigma_sq = ps.dot(entropies**2) - mu**2
    return 2*L - (mu*L + K), sqrt(sigma_sq*L)
    
def validation_plot(L=10, N=50, ref_trials=1000):
    check_points = np.linspace(0,10,10)
    #ics_ref = sorted([motif_ic(random_motif(L, N)) for i in range(ref_trials)])
    ics_ref = sorted([motif_ic(random_motif(L, N), correct=True) for i in trange(ref_trials)])

    plt.plot(ics_ref, 1 - np.linspace(0,1,len(ics_ref)),label="Empirical Complementary CDF", marker='o',linestyle='')
    plt.plot(check_points, [exp(ic_log_pvalue(N, L, ic, method="MC")) for ic in check_points],
             label="Importance Sampling Estimate")
    plt.plot(check_points, [exp(ic_log_pvalue(N, L, ic, method="UB")) for ic in check_points],
             label="Analytic Upper Bound")
    plt.plot(check_points, [exp(ic_log_pvalue(N, L, ic, method="analytic")) for ic in check_points],
             label="Analytic P-value")
    plt.semilogy()
    plt.legend()
    plt.xlabel("Information Content (bits)")
    plt.ylabel("P-value")
    plt.xlim(0,1.2)
    plt.show()
    #return ics_ref

def test_integral(beta, mu, sigma, trials=10000):
    obs = mean(exp(beta*x) for x in [random.gauss(mu, sigma) for i in trange(trials)])
    pred = exp(beta**2*sigma**2/2.0 + beta*mu)
    return pred, obs

def compute_expectation(beta, mu, sigma):
    """compute \int(exp(-beta*x)*phi(x;mu,sigma^2)) from mu to infty"""
    print beta, mu, sigma
    # we do this because the first term is liable to overflow
    second_term = erfc((beta*sigma)/(sqrt(2)))
    if second_term > 0:
        first_term = 1/2.0 * exp(beta**2*sigma**2/2.0 + -beta*mu)
        return first_term * second_term
    else:
        return 0

def compute_expectation_spec(beta, mu, sigma):
    """compute \int(exp(-beta*x)*phi(x;mu,sigma^2)) from mu to infty"""
    print beta, mu, sigma
    # we do this because the first term is liable to overflow
    u = beta*sigma/sqrt(2)
    prefactor = 1/2.0 * exp(-beta*mu)
    # first_term_ref = exp(u**2)
    # second_term_ref = erfc(u)
    # ans_ref = prefactor * first_term * second_term
    ans = prefactor * erfcx(u)
    # print ans_ref, ans
    return ans


def compute_log_expectation(beta, mu, sigma):
    """compute \int(exp(-beta*x)*phi(x;mu,sigma^2)) from mu to infty"""
    print beta, mu, sigma
    # we do this because the first term is liable to overflow
    
    log_second_term = log1p(-erf((beta*sigma**2)/(sigma*sqrt(2))))
    log_first_term = 1/2.0 * exp(beta**2*sigma**2/2.0 + -beta*mu)
    return first_term * second_term

def approx_expectation(beta, mu, sigma, trials=1000):
    return mean(exp(-beta*x) if x > mu else 0 for x in np.random.normal(mu,sigma,trials))
    
def test_expectation():
    beta = random.random() * 10
    mu = random.random() * 10 - 5
    sigma = random.random() * 10
    return compute_expectation(beta, mu, sigma), approx_expectation(beta, mu, sigma, trials=10000)

def alignment_simulation():
    ell = 100
    L = 10
    N = 50
    trials = 100000
    seqs = [random_site(ell) for i in range(N)]
    def random_alignment():
        rs = [random.randrange(ell-L+1) for _ in range(N)]
        return [seq[r:r+L] for seq, r in zip(seqs,rs)]
    ics = [motif_ic(random_alignment()) for _ in trange(trials)]

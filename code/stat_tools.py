import time
import numpy as np
import scipy.stats as stats
from plot_tools import curve_plot

def fisher_test(row1, row2):
    """Calculate statistical significance for a contingency table using Fisher's exact
        test

    Args:
        row1 (list): first row of contingency table
        row2 (list): second row of contingency table
    
    """
    oddsratio, pvalue = stats.fisher_exact([row1, row2])
    print("Fisher's exact test, p-value = %g (odds ratio = %g)" % (pvalue, oddsratio))

def t_test(sample1, sample2):
    """Calculate statistical significance for difference in mean between two data samples
        of unequal variances using two-sided t-test.

    Args:
        sample1 (array): first sample
        sample2 (array): second sample
    
    """
    _, pvalue = stats.ttest_ind(sample1, sample2, equal_var = False)
    print("t-test, p-value = %g" % pvalue)

def perm_test(sample1, sample2, iter):
    """Calculate statistical significance for difference in mean between two data samples
        using permutation test.

    Args:
        sample1 (array): first sample.
        sample2 (array): second sample.
        iter (int): number of resamplings.
    
    """
    n, k = len(sample1), 0
    diff = np.abs(np.mean(sample1) - np.mean(sample2))
    s = np.concatenate([sample1, sample2])
    for j in range(iter):
        np.random.shuffle(s)
        if diff <= np.abs(np.mean(s[:n]) - np.mean(s[n:])):
            k += 1
    print("permutation-test, p-value = %g (%d iterations)" % (k / iter, iter))

def bootstrap_test(sample1, sample2, iter = 10000):
    """Calculate statistical significance for difference in mean between two data samples
        using bootstrap test.

    Args:
        sample1 (array): first sample.
        sample2 (array): second sample.
        iter (int): number of resamplings.
    
    """
    n1, n2, k = len(sample1), len(sample2), 0
    diff = np.abs(np.mean(sample1) - np.mean(sample2))
    sample1 = sample1 - np.mean(sample1)
    sample2 = sample2 - np.mean(sample2)
    for j in range(iter):
        randSample1 = np.random.choice(sample1, size=n1, replace=True)
        randSample2 = np.random.choice(sample2, size=n2, replace=True)
        if diff <= np.abs(np.mean(randSample1) - np.mean(randSample2)):
            k += 1
    print("bootstrap-test, p-value = %g (%d iterations)" % (k / iter, iter))

def normality_test (sample, alpha):
    """Check if a data sample has a normal distribution.

    Args:
        sample (array): sample to be checked.
        alpha (float): cutoff for reporting significant p-value.
    
    """
    _, pvalue = stats.mstats.normaltest(sample)
    print("Normality test, p-value = %g" %pvalue)
    if pvalue > alpha:
        print("Distribution is normal")
    else:
        print("Diftribution is not normal")

def equal_variance_test (sample1, sample2, alpha):
    """Check if two data samples have no statistically significant difference in their
        variances.

    Args:
        sample1 (array): first sample.
        sample2 (array): second sample.
        alpha (float): cutoff for reporting significant p-value.
    
    """
    _, pvalue = stats.levene(sample1, sample2)
    print("Equal variance test, p-value = %g" %pvalue)
    if pvalue < alpha:
        print("Distribution variances are equal")
    else:
        print("Diftribution variances not equal")

def sderror (s):
    """Calculate standard error of the mean.

    Args:
        s (list): population sample.
    
    Returns:
        float: standard error of the mean.
    """
    return np.std(s) / np.sqrt(len(s))

def sderror_on_fraction (k, n):
    """Calculate standard error on a fraction.

    Args:
        k (int): numerator.
        n (int): denominator (sample size).
    
    Returns:
        float: standard error.
    """
    p = k / n
    return np.sqrt(p * (1 - p) / n)

def proportion_ratio_CI (k1, n1, k2, n2):
    
    if ( k1 > 0 ) and ( n1 > 0 ) and ( k2 > 0 ) and ( n2 > 0 ):
        r = (k1 / n1) / (k2 / n2)
        SElogR = np.sqrt( (1 / k1) - (1 / n1) + (1 / k2) - (1 / n2) )
        lowerlog = np.log( r ) - ( 1.96 * SElogR )
        upperlog = np.log( r ) + ( 1.96 * SElogR )
        return np.exp( lowerlog ), np.exp( upperlog )
    else:
        return np.nan, np.nan

def confidence_interval (n,
                         k_obs,
                         pvalues,
                         alpha = 0.05,
                         eventname = 'success',
                         expname = 'experiments',
                         showfig = True,
                         figdir = None,
                         prefigname = 'CI'):
     
    s = time.time()
    prob1, prob2 = [], []
    i1 = i2 = 0
    found1 = found2 = False
    if k_obs > n / 2:
        pvalues = list(reversed(pvalues))
    prev1, prev2 = binomial_prob (n,
                                  k_obs,
                                  pvalues[0])
    for i, p in enumerate(pvalues):
        bp1, bp2 = binomial_prob (n,
                                  k_obs,
                                  p,
                                  itr = 100000)
        if (prev1 <= alpha <= bp1) or (bp1 <= alpha <= prev1):
            i1 = i if abs(bp1 - alpha) <= abs(prev1 - alpha) else i - 1
            found1 = True
        if (prev2 <= alpha <= bp2) or (bp2 <= alpha <= prev2):
            i2 = i if abs(bp2 - alpha) <= abs(prev2 - alpha) else i - 1
            found2 = True
        prob1.append(bp1)
        prob2.append(bp2)
        if found1 and found2:
            break
        prev1 = bp1
        prev2 = bp2
    prob1 = prob1[ : i1 + 1 ]
    prob2 = prob2[ : i2 + 1 ]
    
    for i, prob, (k_min, k_max) in zip( (i1, i2), (prob1, prob2), ((0, k_obs), (k_obs, n))):
        curve_plot ([pvalues[ : i + 1 ]],
                    [prob],
                    styles = '-',
                    xlabel = 'p',
                    ylabel = 'P( %d <= observed %s <= %d )' % (k_min, eventname, k_max),
                    show = showfig,
                    figdir = figdir,
                    figname = '_'.join([prefigname,
                                        'n=%d' % n,
                                        '%d<=k<=%d' % (k_min, k_max)]))
    
    lower_p, upper_p = sorted( ( pvalues[ i1 ], pvalues[ i2 ] )) 
    print('P(≤%d %s events) and P(≥%d %s events) among %d %s ' % (k_obs,
                                                                  eventname,
                                                                  k_obs,
                                                                  eventname,
                                                                  n,
                                                                  expname)
          + 'are both >%.1f%% at %f < p < %f (simulation time = %.1f min)' % (alpha * 100,
                                                                              lower_p,
                                                                              upper_p,
                                                                              (time.time() - s) / 60))
    return lower_p, upper_p

def binomial_prob (n,
                   kobs,
                   p,
                   itr = 10000):
    
    result = []
    k = np.random.binomial(n, p, size = itr)
    for kmin, kmax in [(0, kobs), (kobs, n)]:
        #print('p = %f, %d ≤ k ≤ %d' % (p, kmin, kmax))    
        if p == 0:
            bp = 1 if kmin <= 0 else 0
        elif p == 1:
            bp = 1 if kmax >= n else 0
        else:
            bp = sum( (k >= kmin) & (k <= kmax) ) / itr
        result.append(bp)
    return result

def binomial_test_diff_in_fraction (n1,
                                    n2,
                                    k1,
                                    k2,
                                    itr = 10000):
    """Calculate statistical significance for difference in two fractions by sampling from 
    binomial distributions having equal probabilities of success (p).

    Args:
        n1 (int): number of trials performed on the first distribution.
        n2 (int): number of trials performed on the second distribution.
        k1 (int): number of observed successes from the first distribution.
        k2 (int): number of observed successes from the second distribution.
        itr (int): number of resamplings.
    
    """
    binomial_test_diff_in_fraction_product ([ n1 ],
                                            [ n2 ],
                                            [ k1 ],
                                            [ k2 ],
                                            itr = itr)

def binomial_test_diff_in_fraction_product (n1,
                                            n2,
                                            k1,
                                            k2,
                                            itr = 10000):
    """Calculate statistical significance for difference in two fractions, each being the 
    product of multiple fractions, by sampling from binomial distributions having equal 
    probabilities of success (p).

    Args:
        n1 (list): number of trials performed on each distribution used in the first product.
        n2 (list): number of trials performed on each distribution used in the second product.
        k1 (list): number of observed successes from each distribution used in the first product.
        k2 (list): number of observed successes from each distribution used in the second product.
        itr (int): number of resamplings.
    
    """
    obs_frac1 = np.array( [ k / n for k, n in zip(k1, n1) ] )
    obs_frac2 = np.array( [ k / n for k, n in zip(k2, n2) ] )
    avg_p = ( obs_frac1 + obs_frac2 ) / 2
    obs_diff = np.abs( np.product( obs_frac1 ) - np.product( obs_frac2 ) )
    frac1 = 1
    for f in [ np.random.binomial(n, p, size = itr) / n for n, p in zip(n1, avg_p) ]:
        frac1 *= f
    frac2 = 1
    for f in [ np.random.binomial(n, p, size = itr) / n for n, p in zip(n2, avg_p) ]:
        frac2 *= f
    p_value = sum( np.abs(frac1 - frac2) >= obs_diff ) / itr
    print("binomial sampling test, p-value = %g (%d iterations)" % (p_value, itr))

#----------------------------------------------------------------------------------------
# Modules for statistical computations.
#----------------------------------------------------------------------------------------

import time
import numpy as np
import scipy.stats as stats
from itertools import compress
from collections import Counter

def fisher_test (row1, row2):
    """Calculate statistical significance for a contingency table using Fisher's exact
        test.

    Args:
        row1 (array): first row of contingency table.
        row2 (array): second row of contingency table.
    
    """
    oddsratio, pvalue = stats.fisher_exact([row1, row2])
    print("Fisher's exact test, p-value = %g (odds ratio = %g)" % (pvalue, oddsratio))

def t_test (sample1, sample2):
    """Calculate statistical significance for difference in mean between two data samples
        of unequal variances using two-sided t-test.

    Args:
        sample1 (array): first sample.
        sample2 (array): second sample.
    
    """
    _, pvalue = stats.ttest_ind(sample1, sample2, equal_var = False)
    print("t test, p-value = %g" % pvalue)

def perm_test (sample1, sample2, iter):
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
    print("permutation test, p-value = %g (%d iterations)" % (k / iter, iter))

def bootstrap_test (sample1, sample2, iter = 10000):
    """Calculate statistical significance for difference in mean between two data samples
        by bootstrapping the data samples.

    Args:
        sample1 (array): first sample.
        sample2 (array): second sample.
        iter (int): number of resamplings.
    
    """
    n1, n2, k = len(sample1), len(sample2), 0
    if n1 and n2:
        diff = np.abs(np.mean(sample1) - np.mean(sample2))
        sample1 = sample1 - np.mean(sample1)
        sample2 = sample2 - np.mean(sample2)
        for j in range(iter):
            randSample1 = np.random.choice(sample1, size=n1, replace=True)
            randSample2 = np.random.choice(sample2, size=n2, replace=True)
            if diff <= np.abs(np.mean(randSample1) - np.mean(randSample2)):
                k += 1
        print("bootstrap test, p-value = %g (%d iterations)" % (k / iter, iter))
    else:
        print("bootstrap test, p-value = NaN (empty sample)")

def normality_test (sample, alpha):
    """Check if a data sample has a normal distribution.

    Args:
        sample (array): data sample.
        alpha (float): p-value cutoff for reporting statistical significance.
    
    """
    if len(sample) > 7:
        _, pvalue = stats.mstats.normaltest(sample)
        print("Normality test, p-value = %g" % pvalue)
        if pvalue > alpha:
            print("Distribution is normal")
        else:
            print("Diftribution is not normal")
    else:
        print("Normality test not possible: sample size smaller than 8")

def equal_variance_test (sample1, sample2, alpha):
    """Check if two data samples have no statistically significant difference in their
        variances.

    Args:
        sample1 (array): first sample.
        sample2 (array): second sample.
        alpha (float): p-value cutoff for reporting statistical significance.
    
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
        s (array): data sample.
    
    Returns:
        float

    """
    l = len(s)
    if l:
        return np.std(s) / np.sqrt(l)
    else:
        return np.nan

def sderror_on_fraction (k, n):
    """Calculate standard error on a fraction.

    Args:
        k (int): numerator.
        n (int): denominator (sample size).
    
    Returns:
        float

    """
    p = k / n
    return np.sqrt(p * (1 - p) / n)

def proportion_CI (k, n, conf = 95.0):
    """Calculate the confidence interval on a proportion k/n.

    Args:
        k (int): numerator.
        n (int): denominator (sample size).
        conf (numeric): % confidence interval.
    
    Returns:
        float, float

    """
    mult = { 90.0 : 1.645,
             95.0 : 1.960,
             97.5 : 2.240,
             99.0 : 2.576,
             99.9 : 3.291 }
    if k > 0 and n > 0 and conf in mult:
        se = sderror_on_fraction (k, n)
        p = k / n
        logSE = np.sqrt((1 / k) - (1 / n))
        lowerlog = np.log(p) - mult[conf] * logSE
        upperlog = np.log(p) + mult[conf] * logSE
        return np.exp(lowerlog), np.exp(upperlog)
    else:
        return np.nan, np.nan

def proportion_ratio_CI (k1, n1, k2, n2, a=1, b=0, c=1, d=0, conf = 95.0):
    """Calculate the confidence interval on a ratio of two linear functions of proportions 
        (a * (k1/n1) + b) / (c * (k2/n2) + d).

    Args:
        k1, k2, n1, n2 (int): proportion numerators and denominators.
        a, b, c, d (float): prportion linear function parameters.
        conf (numeric): % confidence interval.
    
    Returns:
        float, float

    """
    mult = { 90.0 : 1.645,
             95.0 : 1.960,
             97.5 : 2.240,
             99.0 : 2.576,
             99.9 : 3.291 }
    if k1 > 0 and n1 > 0 and k2 > 0 and n2 > 0 and conf in mult:
        r = (a * k1/n1 + b) / (c * k2/n2 + d)
        var1 = (1/k1 - 1/n1) /  (1 + (b * n1) / (a * k1))
        var2 = (1/k2 - 1/n2) / (1 + (d * n2) / (c * k2))
        SElogR = np.sqrt(var1 + var2)
        lowerlog = np.log(r) - mult[conf] * SElogR
        upperlog = np.log(r) + mult[conf] * SElogR
        return np.exp(lowerlog), np.exp(upperlog)
    else:
        return np.nan, np.nan

def proportion_sum_CI (k1, n1, k2, n2, a=1, b=1, conf = 95.0):
    """Calculate the confidence interval on a sum of two proportions 
        a(k1/n1) + b(k2/n2).

    Args:
        k1, k2, n1, n2 (int): proportion numerators and denominators.
        a, b (float): prportion coefficients.
        conf (numeric): % confidence interval.
    
    Returns:
        float, float

    """
    mult = { 90.0 : 1.645,
             95.0 : 1.960,
             97.5 : 2.240,
             99.0 : 2.576,
             99.9 : 3.291 }
    if k1 > 0 and n1 > 0 and k2 > 0 and n2 > 0 and conf in mult:
        p1 = k1 / n1
        p2 = k2 / n2
        s = a*p1 + b*p2
        term1 = a**2 * p1 * (1 - p1) / n1
        term2 = b**2 * p2 * (1 - p2) / n2
        term3 = a*p1 + b*p2
        SElog = np.sqrt(term1 + term2) / term3
        lowerlog = np.log(s) - mult[conf] * SElog
        upperlog = np.log(s) + mult[conf] * SElog
        return np.exp(lowerlog), np.exp(upperlog)
    else:
        return np.nan, np.nan

def confidence_interval (n,
                         k_obs,
                         pvalues,
                         alpha = 0.05,
                         eventname = 'success',
                         expname = 'experiments'):
    """Calculate the confidence interval on a proportion k_obs/n using a Binomial 
        distribution produced from Bernoulli trials.

    Args:
        n (int): proportion denominator.
        k_obs (int): proportion numerator.
        pvalues (array): probabilities of success used for Bernoulli trials.
        alpha (float): p-value cutoff for reporting statistical significance.
        eventname (str): name of success event in the Bernoulli trial (for display purpose).
        expname (str): name of experiment on which Bernoulli trials are performed (for display purpose).
    
    Returns:
        float, float

    """
    s = time.time()
    prob1, prob2 = [], []
    i1 = i2 = 0
    found1 = found2 = False
    if k_obs > n / 2:
        pvalues = list(reversed(pvalues))
    prev1, prev2 = binomial_prob (n, k_obs, pvalues[0])
    for i, p in enumerate(pvalues):
        bp1, bp2 = binomial_prob (n, k_obs, p, itr = 100000)
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
        prev1, prev2 = bp1, bp2
    prob1 = prob1[:i1 + 1]
    prob2 = prob2[:i2 + 1]
    
    lower_p, upper_p = sorted( ( pvalues[i1], pvalues[i2] )) 
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

def binomial_prob (n, kobs, p, itr = 10000):
    """Calculate the probability of observing ≤ kobs or ≥ kobs successes in a Bernoulli 
        experiment with probability of success p.

    Args:
        n (int): number of trials in experiment.
        kobs (int): number of successes observed in experiment.
        p (float): probability of success used in Bernoulli trials.
        itr (int): number of Bernoulli trials.

    Returns:
        float, float

    """
    result = []
    k = np.random.binomial(n, p, size = itr)
    for kmin, kmax in [(0, kobs), (kobs, n)]:
        #print('p = %f, %d ≤ k ≤ %d' % (p, kmin, kmax))    
        if p == 0:
            bp = 1 if kmin <= 0 else 0
        elif p == 1:
            bp = 1 if kmax >= n else 0
        else:
            bp = sum((k >= kmin) & (k <= kmax)) / itr
        result.append(bp)
    return result

def binomial_test_diff_in_fraction (n1, n2, k1, k2, itr = 10000):
    """Calculate statistical significance for difference in two fractions by sampling from 
    binomial distributions having equal probabilities of success (p).

    Args:
        n1 (int): number of trials performed on the first distribution.
        n2 (int): number of trials performed on the second distribution.
        k1 (int): number of observed successes from the first distribution.
        k2 (int): number of observed successes from the second distribution.
        itr (int): number of resamplings.
    
    """
    binomial_test_diff_in_fraction_product ([n1], [n2], [k1], [k2], itr = itr)

def binomial_test_diff_in_fraction_product (n1, n2, k1, k2, itr = 10000):
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

def compress_data (xpos, ypos, xstart = 0, binwidth = 0.1, perbin = 0):
    """Compress y-axis data based on x-axis binning.

    Args:
        xpos (array): x-axis data points.
        ypos (array): y-axis data points.
        xstart (numeric): starting point of first bin on x-axis.
        binwidth (numeric): bin width on x-axis.
        perbin (int): maximum number of data points compressed per bin.

    Returns:
        float, float, float, float, list: x-data means, y-data means, x-data SEM, y-data SEM, bin borders.

    """
    xposmax = max(xpos)
    xmin, xmax = xstart, xstart + binwidth
    xmean, ymean, xerr, yerr, bins = [], [], [], [], []
    while xmin < xposmax:
        if xmin == xstart:
            select = [ xmin <= x <= xmax for x in xpos ]
        else:
            select = [ xmin < x <= xmax for x in xpos ]
        numpts = sum(select)
        if (numpts == 0) and (perbin == 0):
            xmin = xmax
            xmax = xmin + binwidth
        elif (numpts >= perbin) or (xmax >= xposmax):
            xcompress = list( compress(xpos, select) )
            ycompress = list( compress(ypos, select) )
            bins.append( (xmin, xmax) )
            xmean.append( np.mean( xcompress ) )
            ymean.append( np.mean( ycompress ) )
            xerr.append( sderror( xcompress ) )
            yerr.append( sderror( ycompress ) )
            xmin = xmax
            xmax = xmin + binwidth
        else:
            xmax = xmax + binwidth
    return xmean, ymean, xerr, yerr, bins

def round_data (data, w = 1, maxVal = np.inf):
    """Round data points.

    Args:
        data (array): data points to be rounded.
        w (numeric): rounding interval.
        maxVal (numeric): maximum rounded value.

    Returns:
        list: rounded data.

    """
    rounded = []
    for d in data:
        if d <= maxVal - w/2:
            rounded.append(w * round(float(d) / w))
        else:
            rounded.append(maxVal)
    return rounded

def pdf (data):
    """Produce a discrete density distribution for a sample of data.

    Args:
        data (array): data points.

    Returns:
        array, array: density distribution, unique data points.

    """
    dataCount = Counter(data)
    x = np.array(sorted(dataCounter.keys()))
    dens = np.array( [dataCount[d] for d in x] )    
    return dens/sum(dens), x

def cont_pdf (data, minVal = None, maxVal = None, w = 1):
    """Produce a continuous density distribution for a sample of data.

    Args:
        data (array): data points.
        minVal (numeric): starting point of density distribution.
        maxVal (numeric): ending point of density distribution.
        w (numeric): interval between data point.

    Returns:
        array, array: density distribution, unique data points.

    """
    dataCount = Counter(data)
    k = dataCount.keys()
    if minVal is None:
        minVal = min(k)
    if maxVal is None:
        maxVal = max(k)
    x = np.arange(minVal, maxVal + w, w)
    dens = np.array( [dataCount[d] for d in x] )    
    return dens/sum(dens), x

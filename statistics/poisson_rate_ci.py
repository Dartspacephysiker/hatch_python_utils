import scipy.stats as stats

def prci(N,T,
                    alpha=0.05):
    """
    rate_L,rate_U = prci(N,alpha=0.05)
    
    INPUT
    =====
    N            (int)      : Number of observations
    T            (int/float): Total amount of time for which observations were made

    KEYWORDS
    ========
    alpha        (float)    : 1-[desired confidence interval] (e.g., alpha=0.05 corresponds to the 95% confidence interval)
    

    OUTPUT
    ======
    rate_L                  : Lower bound of the rate for the selected confidence interval
    rate_U                  : Upper bound of the rate for the selected confidence interval

    REFERENCES
    ==========
    •Ulm K. A simple method to calculate the confidence interval of a standardized mortality ratio (SMR). Am J Epidemiol. 1990 Feb;131(2):373-5. doi: 10.1093/oxfordjournals.aje.a115507. PMID: 2296988.
    •https://www.statsdirect.com/help/rates/poisson_rate_ci.htm
    
    AUTHOR
    ======
    Spencer Mark Hatch
    University of Bergen
    2020-09-24

    EXAMPLE (from https://www.statsdirect.com/help/rates/poisson_rate_ci.htm):
    =======

    # Inputs from website
    N,T = 14,400

    # 95% CI from website
    rate_L_CORRECT,rate_U_CORRECT = 0.019135,0.058724

    # Calculate 95% CI
    rate_L,rate_U = prci(N,T,alpha=0.05)

    print("rate_L          : {:8.5f}".format(rate_L))
    print("rate_L (correct): {:8.5f}".format(rate_L_CORRECT))
    print("")
    print("rate_U          : {:8.5f}".format(rate_U))
    print("rate_U (correct): {:8.5f}".format(rate_U_CORRECT))
    """

    rate = N/T

    dfL = N*2
    dfU = (N+1)*2

    # Use "percent point function" to get these values
    Nl = stats.chi2.ppf(alpha/2,dfL)/2
    Nu = stats.chi2.ppf(1-alpha/2,dfU)/2

    rateL,rateU = Nl/T,Nu/T

    return rateL,rateU


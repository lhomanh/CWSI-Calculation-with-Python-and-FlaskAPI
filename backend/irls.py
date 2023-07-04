from numba import jit_module
from numpy import abs, diag, ones, sort, where
from numpy.linalg import inv

CACHE = False

def median(a):
    '''median value


    Arguments
    ---------
    a : float
        array of values


    Returns
    -------
    float
    '''
    a = sort(a)
    n = len(a)
    m = int(n / 2) # floor division
    if n % 2: # if odd
        return a[m]
    else:
        return (a[m] + a[m - 1]) / 2

def design(x, deg=1):
    '''polynomial design matrix


    Arguments
    ---------
    x : float (n)
        independent values
    deg : int, optional
        polynomial degree


    Returns
    -------
    float (n, deg + 1)
    '''
    n = len(x)
    res = ones((n, deg + 1))
    for i in range(1, deg + 1):
        res[:,i] = x ** i
    return res

def yhat(b, x):
    '''estimated y


    Arguments
    ---------
    b : float (m)
        regression coefficients
    x : float (n)
        independent values


    Returns
    -------
    float (n)
    '''
    x = design(x, len(b) - 1)
    return x @ b

def residuals(X, b, y):
    '''robust residuals


    Arguments
    ---------
    X : float (n, m)
        independent design matrix
    b : float (m)
        regression coefficients
    y : float (n)
        dependent values


    Returns
    -------
    float (n)
    '''
    n = len(y)
    e = y - X @ b # errors
    # median of absolute deviations from the median error
    mad = median(abs(e - median(e))) * 0.6745
    return e / mad # scaled residuals

def wls(X, w, y):
    '''weighted least squares regression coefficients


    Arguments
    ---------
    X : float (n, m)
        independent design matrix
    w : float (n)
        regression coefficients
    y : float (n)
        dependent values


    Returns
    -------
    float (m)
    '''
    Xtw = X.T @ diag(w) # weighted X transpose
    return inv(Xtw @ X) @ Xtw @ y # regression parameters

def tukey(u):
    '''Tukey's bisquare weights


    Arguments
    ---------
    u : float (n)
        residuals

    Returns
    -------
    float (n)
    '''
    return where(abs(u) > 4.685, 0, (1 - (u / 4.685) ** 2) ** 2)

def huber(u):
    '''Huber weights


    Arguments
    ---------
    u : float (n)
        residuals

    Returns
    -------
    float (n)
    '''
    u = abs(u)
    return where(u > 1.345, 1.345 / u, 1)

def irls(X, y, limit=666, th=0.001, weight=tukey):
    '''iteratively reweighted lease squares regression coefficients with final
    weights


    Arguments
    ---------
    X : float (n, m)
        independent design matrix
    y : float (n)
        dependent values
    limit : int, optional
        maximum iterations
    th : float, optional
        relative convergence threshold for all coefficients
    weight : function, optional
        weight function, None returns OLS solution


    Returns
    -------
    float (m)
        coefficients
    float (n)
        final weights
    '''
    w = ones(len(X)) # initial weight array
    b = wls(X, w, y) # initial OLS regression
    if not weight:
        return b, w
    u = residuals(X, b, y)
    w = weight(u) # initial weights
    itr = 0 # count interations
    while itr < limit: # limit iterations to limit
        itr += 1 # increment counter
        tmp = wls(X, w, y) # first weighted run
        if (abs(b - tmp) / b <= th).all(): # if all parameters in tolerance
            return tmp, w # return parameters
        b = tmp # else tmp is new set of parameters
        u = residuals(X, b, y)
        w = weight(u) # get new weights and repeat
    assert False

#jit_module(nopython=True, inline='always', cache=CACHE)

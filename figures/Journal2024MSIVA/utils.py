import scipy
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.linear_model import Ridge
import matplotlib.pyplot as plt

def zscore(data, axis):
    data -= data.mean(axis=axis, keepdims=True)
    data /= data.std(axis=axis, keepdims=True)
    return np.nan_to_num(data, copy=False)


def correlation(matrix1, matrix2):
    d1 = matrix1.shape[-1]
    d2 = matrix2.shape[-1]

    assert d1 == d2
    assert matrix1.ndim <= 2
    assert matrix2.ndim <= 2
    
    matrix1 = zscore(matrix1.astype(float), matrix1.ndim - 1) / np.sqrt(d1)
    matrix2 = zscore(matrix2.astype(float), matrix2.ndim - 1) / np.sqrt(d2)
    
    if matrix1.ndim >= matrix2.ndim:
        return np.dot(matrix1, matrix2.T)
    else:
        return np.dot(matrix2, matrix1.T)


def compute_rdc(X, Y):
    d = X.shape[0]
    cc = np.zeros((d, d))
    for i in range(d):
        for j in range(d):
            cc[i, j] = np.abs(rdc(X[i, :], Y[j, :]))
    return cc


def rdc(x, y, k=10, s=.5, nonlinearity='sin', seed=0):
    """
    Python implementation of the Randomized Dependence Coefficient (RDC) [1] algorithm
    the RDC is a measure of correlation between two (scalar) random variables x and y
    that is invariant to permutation, scaling, and most importantly nonlinear scaling

    Parameters:
        x: numpy array of shape (n,)
        y: numpy array of shape (n,)
        k: number of random projections in RDC
        s: covariance of the Gaussian dist used for sampling the random weights
        nonlinearity: nonlinear feature map used to transform the random projections

    Return:
        rdc_cc: flaot in [0,1] --- the RDC correlation coefficient

    References:
    [1] https://papers.nips.cc/paper/2013/file/aab3238922bcc25a6f606eb525ffdc56-Paper.pdf

    """
    cx = copula_projection(x, k, s, nonlinearity, seed)
    cy = copula_projection(y, k, s, nonlinearity, seed)
    rdc_cc = largest_cancorr(cx, cy)
    return rdc_cc


def copula_projection(x, k=20, s=.5, nonlinearity='sin', seed=0):
    n = x.shape[0]
    k = min(k, n)
    # compute the empirical cdf (copula) of x evaluated at x
    p = rank_array(x) / n  # (n, )
    # augment the copula with 1
    pt = np.vstack([p, np.ones(n)]).T  # (n, 2)
    # sample k random weights
    np.random.seed(seed)
    wt = np.random.normal(0, s, size=(pt.shape[1], k))
    if nonlinearity == 'sin':
        phix = np.sin(pt.dot(wt))  # (n, k)
    elif nonlinearity == 'cos':
        phix = np.cos(pt.dot(wt))  # (n, k)
    else:
        raise ValueError(f'{nonlinearity} not supported')
    return np.hstack([phix, np.ones((n, 1))])


def make_diag(el, nrows, ncols):
    diag = np.zeros((nrows, ncols))
    for i in range(min(nrows, ncols)):
        diag[i, i] = el
    return diag


def largest_cancorr(x, y):
    """
    Return the largest correlation coefficient after solving CCA between two matrices `x` and `y`.
    inspired from R's `cancor` function
    """
    n = x.shape[0]
    x = x - x.mean(axis=0)
    y = y - y.mean(axis=0)
    qx, _ = scipy.linalg.qr(x, mode='full')
    qy, _ = scipy.linalg.qr(y, mode='full')
    dx = np.linalg.matrix_rank(x)
    dy = np.linalg.matrix_rank(y)
    qxy = qx.T.dot(qy.dot(make_diag(1, n, dy)))[:dx]
    _, s, _ = scipy.linalg.svd(qxy, lapack_driver='gesvd')
    return s[0]


def rank_array(x):
    tmp = x.argsort()
    ranks = np.empty_like(tmp)
    ranks[tmp] = np.arange(len(x))
    return ranks + 1


def age_regression(X_train, X_test, y_train, y_test, a = 1):
    clf = Ridge(alpha=a) 
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    mae = np.mean(np.abs(y_pred-y_test))
    return mae


def calculate_balanced_accuracy(y_true, y_pred):
    n_classes = len(np.unique(y_true))
    acc = np.zeros(n_classes)
    for i in range(n_classes):
        acc[i] = np.mean(y_pred[y_true == i] == i)
    return np.mean(acc)


def sex_classification(X_train, X_test, y_train, y_test, c=1):
    clf = LinearSVC(C=c,dual=False)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    # acc = 1 - np.mean(np.abs(y_pred-y_test))
    acc = calculate_balanced_accuracy(y_test, y_pred)
    return acc


def sz_classification(X_train, X_test, y_train, y_test, c=1):
    clf = LinearSVC(C=c,dual=False)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    # acc = 1 - np.mean(np.abs(y_pred-y_test))
    acc = calculate_balanced_accuracy(y_test, y_pred)
    return acc


def plot_sq(ax, i, crossmodal=False):
    sqc1 = "palegreen"
    sqc2 = "cornflowerblue"
    s_dict = {0:{0:2, 2:3, 5:4, 9:1, 10:1, 11:1},\
            1:{0:2, 2:2, 4:2, 6:2, 8:2, 10:1, 11:1},\
            2:{0:3, 3:3, 6:3, 9:1, 10:1, 11:1},\
            3:{0:4, 4:4, 8:1, 9:1, 10:1, 11:1},\
            4:{0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1, 10:1, 11:1}}
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sq_dict = s_dict.get(i)
    for k in sq_dict.keys():
        v = sq_dict[k]
        if crossmodal and v == 1 and i != 4:
            continue
        if v == 1 and i != 4:
            sqc=sqc1
        else:
            sqc=sqc2
        square = plt.Rectangle((k,k), v, v, fill=False, ec=sqc, linewidth=2, zorder=10)
        ax.add_patch(square)
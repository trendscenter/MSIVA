import scipy
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.svm import LinearSVC
from sklearn.linear_model import Ridge


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


def compute_itc(X, n_source):
    """
    Compute information criteria from multimodal group PCA.
    :param X: list of numpy arrays, each of shape (n_features, n_samples)
    :param n_source: int, number of components to extract
    :return: aic: numpy array, shape (n_source-1,)
    :return: kic: numpy array, shape (n_source-1,)
    :return: mdl: numpy array, shape (n_source-1,)
    """
    M = range(len(X))
    n_sample = X[0].shape[1]
    V = [x.shape[0] for x in X]
    
    if n_sample <= min(V):
        cvx = np.zeros((n_sample, n_sample))
        
        for mm in M:
            cvx_ = np.cov(X[mm], rowvar=False)
            cvx += cvx_ / (len(M) * np.trace(cvx_) / n_sample)
        
        # Subject-level PCA reduction...
        lambda_, H = np.linalg.eig(cvx)
        lambda_ = lambda_[:n_source]
        H = H[:, :n_source]
    
    else:
        X_cat = np.concatenate([np.sqrt(n_sample / (len(M) * np.sum(x**2))) * x for x in X], axis=0)
        cvx = X_cat @ X_cat.T
        
        lambda_, U = np.linalg.eig(cvx)
        lambda_ = lambda_[:n_source]
        U = U[:, :n_source]
        H = ((np.diag(1.0 / np.sqrt(lambda_)) @ U.T) @ X_cat).T
    
    # Step 1: Compute U and normalize
    U = [(x @ H) for x in X]
    norm_U = [np.sum(u**2, axis=0) for u in U]
    norm_U = np.sqrt(np.sum(np.array(norm_U), axis=0))
    U = [u / norm_U for u in U]
    U = np.concatenate(U, axis=1).T

    # Step 2: Compute square root of eigenvalues
    lambda_sqrt = np.sqrt(lambda_)

    aic = np.zeros(n_source - 1)
    kic = np.zeros(n_source - 1)
    mdl = np.zeros(n_source - 1)

    # Step 3: Loop through and compute AIC, KIC, and MDL
    for k in range(n_source - 1):
        LH = (1 / (n_source - k)) * np.sum(np.log(lambda_sqrt[k+1:])) - np.log(np.mean(lambda_sqrt[k+1:]))
        mlh = 0.5 * n_sample * (n_source - k) * LH
        df = 1 + 0.5 * k * (2 * n_source - k + 1)
        
        aic[k] = -2 * mlh + 2 * df
        kic[k] = -2 * mlh + 3 * df
        mdl[k] = -mlh + 0.5 * df * np.log(n_sample)
    
    return aic, kic, mdl


def calculate_rdc(X, Y):
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
        rdc_cc: float in [0,1] --- the RDC correlation coefficient

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


def age_regression(X_train, X_test, y_train, y_test, a=1):
    clf = Ridge(alpha=a) 
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    mae = np.mean(np.abs(y_pred-y_test))
    return mae, clf.coef_


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
    return acc, clf.coef_


def sz_classification(X_train, X_test, y_train, y_test, c=1):
    clf = LinearSVC(C=c,dual=False)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    # acc = 1 - np.mean(np.abs(y_pred-y_test))
    acc = calculate_balanced_accuracy(y_test, y_pred)
    return acc, clf.coef_


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
    for k, v in sq_dict.items():
        if crossmodal and v == 1 and i != 4:
            continue
        if v == 1 and i != 4:
            sqc=sqc1
        else:
            sqc=sqc2
        square = plt.Rectangle((k, k), v, v, fill=False, ec=sqc, linewidth=2, zorder=10)
        ax.add_patch(square)


def myISI(WAr):
    N = WAr.shape[0]
    WAr = np.abs(WAr)
    ISI = 0.
    ISI += np.sum(np.sum(WAr,axis=1)/np.max(WAr,axis=1) - 1)
    # np.max(WAr,axis=1)
    ISI += np.sum(np.sum(WAr,axis=0)/np.max(WAr,axis=0) - 1)
    # np.max(WAr,axis=0)
    ISI = ISI/(2*N*(N-1))
    return ISI


def calculate_misi(W,A,S):
    WA = [W[cc] @ A[cc] for cc in range(len(W))]
    
    M = len(S)
    K = S[0].shape[0]

    WAr = np.zeros((K,K))

    for mm in range(M):
        mk = np.split(S[mm] == 1, K, axis=0)
        for kr in range(K):        # rows
            mkr = mk[kr]
            for kc in range(K):    # columns
                mkc = mk[kc]
                WAr[kr,kc] += np.abs(WA[mm][mkr.T @ mkc]).sum()
  
    out = myISI(WAr)

    return out, WAr


def find_all_blocks(matrix, block_sizes):
    """
    Find all possible blocks of different sizes in the matrix.
    Returns a list of (block, row_indices, col_indices, size) tuples.
    """
    n = matrix.shape[0]
    all_blocks = []
    
    for size in block_sizes:
        for i in range(n - size + 1):
            for j in range(n - size + 1):
                block = matrix[i:i+size, j:j+size]
                row_indices = list(range(i, i+size))
                col_indices = list(range(j, j+size))
                all_blocks.append((block, row_indices, col_indices, size))
    
    return all_blocks


def calculate_block_coherence(block):
    """
    Calculate a coherence score for a block.
    Higher scores indicate more structured/coherent blocks.
    """
    # Multiple criteria for block coherence
    variance = np.var(block)
    mean_abs = np.mean(np.abs(block))
    max_val = np.max(np.abs(block))
    
    # Combine criteria (you can adjust these weights)
    coherence_score = variance + mean_abs + max_val
    return coherence_score


def find_non_overlapping_blocks(all_blocks, target_sizes=[2, 3, 4]):
    """
    Find the best non-overlapping combination of blocks with target sizes.
    """
    # Group blocks by size
    blocks_by_size = {}
    for block_data in all_blocks:
        size = block_data[3]
        if size not in blocks_by_size:
            blocks_by_size[size] = []
        blocks_by_size[size].append(block_data)
    
    # Sort blocks by coherence score for each size
    for size in blocks_by_size:
        blocks_by_size[size].sort(key=lambda x: calculate_block_coherence(x[0]), reverse=True)
    
    best_combination = None
    best_score = -1
    
    # Try different combinations of blocks
    for block_2 in blocks_by_size.get(2, []):
        for block_3 in blocks_by_size.get(3, []):
            for block_4 in blocks_by_size.get(4, []):
                if not blocks_overlap([block_2, block_3, block_4]):
                    total_score = (calculate_block_coherence(block_2[0]) +
                                 calculate_block_coherence(block_3[0]) +
                                 calculate_block_coherence(block_4[0]))
                    if total_score > best_score:
                        best_score = total_score
                        best_combination = [block_2, block_3, block_4]
    
    return best_combination if best_combination else []


def blocks_overlap(blocks):
    """
    Check if any blocks overlap in their row or column indices.
    """
    used_rows = set()
    used_cols = set()
    
    for _, row_indices, col_indices, _ in blocks:
        for row in row_indices:
            if row in used_rows:
                return True
            used_rows.add(row)
        for col in col_indices:
            if col in used_cols:
                return True
            used_cols.add(col)
    
    return False


def create_diagonal_permutation(blocks, matrix_size=9):
    """
    Create row and column permutations to place blocks on diagonal.
    Blocks are arranged in order: 2x2, 3x3, 4x4
    """
    # Sort blocks by size for consistent diagonal placement
    blocks = sorted(blocks, key=lambda x: x[3])
    
    row_mapping = {}  # old_row -> new_row
    col_mapping = {}  # old_col -> new_col
    
    current_diagonal_pos = 0
    
    for block, row_indices, col_indices, size in blocks:
        # Map original block position to diagonal position
        for i, old_row in enumerate(row_indices):
            new_row = current_diagonal_pos + i
            row_mapping[old_row] = new_row
            
        for i, old_col in enumerate(col_indices):
            new_col = current_diagonal_pos + i
            col_mapping[old_col] = new_col
        
        current_diagonal_pos += size
    
    # Handle remaining rows and columns
    used_rows = set(row_mapping.keys())
    used_cols = set(col_mapping.keys())
    used_new_rows = set(row_mapping.values())
    used_new_cols = set(col_mapping.values())
    
    remaining_old_rows = [i for i in range(matrix_size) if i not in used_rows]
    remaining_old_cols = [i for i in range(matrix_size) if i not in used_cols]
    remaining_new_rows = [i for i in range(matrix_size) if i not in used_new_rows]
    remaining_new_cols = [i for i in range(matrix_size) if i not in used_new_cols]
    
    # Map remaining rows and columns
    for old_row, new_row in zip(remaining_old_rows, remaining_new_rows):
        row_mapping[old_row] = new_row
    for old_col, new_col in zip(remaining_old_cols, remaining_new_cols):
        col_mapping[old_col] = new_col
    
    return row_mapping, col_mapping


def apply_permutation_from_mapping(matrix, row_mapping, col_mapping):
    """
    Apply row and column permutations using mapping dictionaries.
    """
    n = matrix.shape[0]
    
    # Create permutation arrays
    row_order = [0] * n
    col_order = [0] * n
    
    for old_pos, new_pos in row_mapping.items():
        row_order[new_pos] = old_pos
    
    for old_pos, new_pos in col_mapping.items():
        col_order[new_pos] = old_pos
    
    # Apply permutations
    permuted_matrix = matrix[row_order][:, col_order]
    
    return permuted_matrix, row_order, col_order


def sort_mixed_blocks_matrix(matrix):
    """
    Main function to sort 9x9 matrix with 2x2, 3x3, and 4x4 blocks on diagonal.
    """
    block_sizes = [2, 3, 4]
    
    # Find all possible blocks
    all_blocks = find_all_blocks(matrix, block_sizes)
    # print(f"Found {len(all_blocks)} total possible blocks")
    
    # Find the best non-overlapping combination
    selected_blocks = find_non_overlapping_blocks(all_blocks)
    
    if len(selected_blocks) != 3:
        print(f"Warning: Could only find {len(selected_blocks)} non-overlapping blocks")
        return matrix, list(range(9)), list(range(9)), selected_blocks
    
    # print("Selected blocks:")
    # for i, (block, row_idx, col_idx, size) in enumerate(selected_blocks):
    #     print(f"  {size}x{size} block at rows {row_idx}, cols {col_idx}")
    #     print(f"  Coherence score: {calculate_block_coherence(block):.3f}")
    
    # Create permutations
    row_mapping, col_mapping = create_diagonal_permutation(selected_blocks)
    
    # Apply permutations
    sorted_matrix, row_order, col_order = apply_permutation_from_mapping(
        matrix, row_mapping, col_mapping)
    
    return sorted_matrix, row_order, col_order, selected_blocks


def calculate_mcc(corr, ss, sort=True, col_order=None):
    total_subspace_size = sum(ss)
    abscorr = np.abs(corr[:total_subspace_size, :total_subspace_size])
    
    if sort:
        # print("Sorting the correlation matrix according to subspace structure...")
        if ss[0] == 2 and ss[1] == 3 and ss[2] == 4: # subspace structures with a 2x2 block, a 3x3 block, and a 4x4 block
            sorted_corr, row_order, col_order, found_blocks = sort_mixed_blocks_matrix(abscorr)
        elif col_order is None: # subspace structures with same block size
            K = len(ss)
            aggcorr = np.zeros((K,K))
            for k1, d1 in enumerate(ss):
                ind_start1 = sum(ss[:k1])
                for k2, d2 in enumerate(ss):
                    ind_start2 = sum(ss[:k2])
                    abscorr_subspace = abscorr[ind_start1:ind_start1+d1,ind_start2:ind_start2+d2]
                    aggcorr[k1, k2] = np.mean(abscorr_subspace)
            ind = np.argmax(np.abs(aggcorr),axis=1)
            col_order = np.zeros(np.sum(ss))
            block_size = ss[0]
            for i, j in enumerate(ind):
                col_order[i*block_size:(i+1)*block_size] = np.arange(j*block_size, (j+1)*block_size)
            col_order = col_order.astype(int)
            sorted_corr = abscorr[:,col_order]
        elif col_order is not None:
            sorted_corr = abscorr[:,col_order]
    else:
        sorted_corr = abscorr

    K = len(ss)
    aggcorr = np.zeros((K,K))
    for k1, d1 in enumerate(ss):
        ind_start1 = sum(ss[:k1])
        for k2, d2 in enumerate(ss):
            ind_start2 = sum(ss[:k2])
            abscorr_subspace = sorted_corr[ind_start1:ind_start1+d1,ind_start2:ind_start2+d2]
            maxrow = np.max(abscorr_subspace, axis=1)
            maxcol = np.max(abscorr_subspace, axis=0)
            aggcorr[k1, k2] = (np.sum(maxrow) + np.sum(maxcol))/(d1+d2)
            
    mcc = np.mean(np.diag(aggcorr))
    off_diag_ind = ~np.eye(aggcorr.shape[0], dtype=bool)
    md = (np.sum(aggcorr[off_diag_ind]) + np.sum(1-np.diag(aggcorr)))/(K**2)
    # print(f"mcc: {mcc}, md: {md}")
    return mcc, md, aggcorr, col_order


def min_max_aggregate(aggcorr1, aggcorr2):
    aggcorr = np.zeros_like(aggcorr1)
    for i in range(aggcorr1.shape[0]):
        for j in range(aggcorr1.shape[1]):
            if i == j:
                aggcorr[i,j] = np.min([aggcorr1[i,j], aggcorr2[i,j]])
            else:
                aggcorr[i,j] = np.max([aggcorr1[i,j], aggcorr2[i,j]])
    return aggcorr


def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"


def convert_pvalue(p):
    if p < 0.001:
        return f"<0.001"
    elif p < 0.01:
        return f"<0.010"
    elif p < 0.05:
        return f"<0.050"
    elif p > 0.05:
        return f">0.050"


def pvalue_to_r(p_value, n):
    t_stat = stats.t.ppf(1 - p_value / 2, n - 2)
    r = t_stat / np.sqrt(t_stat**2 + (n - 2))
    return r
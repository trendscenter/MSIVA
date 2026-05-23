#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 2 20:19:00 2018 (MDT)
@author: Rogers F. Silva
Copied from: https://github.com/rsilva8/decentralized-pca/blob/master/ancillary.py
"""

import numpy as np
from scipy import linalg


def base_PCA(data, num_PC=None, axis=-2, whitening=True):
    """ Efficient principal component analysis of a dataset
    The transformation (and reduction) occurs along the selected axis.
    Does not perform data mean removal prior to PCA.
    Does not perform data scaling prior to PCA.
    Parameters
    ----------
    data : 2D array
        The data to be transformed and decomposed.
    num_PC : positive int, optional
        The number of principal components to be retained. 'None' is interpreted as unspecified and
        yields all principal components (up to the matrix rank of data). The default is None.
    axis : negative int, optional
        Axis along which the transformation (and reduction) is applied. The default is -2.
    whitening : bool, optional
        If this is set to True, the principal components are standardized, i.e., the
        principal components are rescaled to have unit standard deviation.
    Returns
    -------
    reduced_data : 2D array
        The transformed (and reduced) dataset (aka, the principal components).
    projM : 2D array
        The matrix that projects the data onto the subspace (aka, the whitening matrix)
    bkprojM : 2D array
        The matrix that reconstructs the data from the transformed (and reduced) dataset
        (aka, the dewhitening or back-projection matrix).
    Notes
    -----
    Implemented based on the eigenvalue decomposition (EVD) of the non-normalized covariance matrix.
    The covariance and EVD are always computed on the shortest axis of the input ``data``.
    The incomplete EVD is performed using ``scipy.linalg.eigh()`` for efficiency,
    unless the covariance matrix is ill-conditioned.
    """

    # base_PCA() must be available/defined both at the remote and local sites.

    num_rows, num_cols = data.shape

    # Compute covariance matrix C along smallest dimension
    if num_rows <= num_cols:
        C = data @ data.T
    else:
        C = data.T @ data

    U, S = do_cov_EVD(C, k=num_PC)

    del C

    # Reduce selected dimension/axis
    if num_rows <= num_cols:
        U = U.T
        tmp = U @ data
        if axis == -2:
            reduced_data = tmp
            projM = U
        elif axis == -1:
            tmp = tmp.T * (1 / np.sqrt(S))
            reduced_data = data @ tmp
            projM = tmp
    else:
        tmp = data @ U
        if axis == -2:
            tmp = (1 / np.sqrt(S[:, None])) * tmp.T
            reduced_data = tmp @ data
            projM = tmp
        elif axis == -1:
            reduced_data = tmp
            projM = U

    del tmp
    del U, S

    # Apply whitening if required.
    # Whitening occurs along the "other" dimension, not the reduced dimension
    # Whitening changes the std dev to 1 by multiplying by the inverse of its std dev
    if whitening:
        if axis == -2:
            stds = np.std(reduced_data, axis=-1, ddof=1)
            invs = 1 / stds
            return invs[:,
                        None] * reduced_data, invs[:,
                                                   None] * projM, projM.T * stds
        elif axis == -1:
            stds = np.std(reduced_data, axis=-2, ddof=1)
            invs = 1 / stds
            return reduced_data * invs, projM * invs, stds[:, None] * projM.T
    else:
        return reduced_data, projM, projM.T


def do_cov_EVD(C, k=None, method=0):
    ''' Compute top ``k`` largest eigenvectors of (symmetric and real-valued) covariance ``C``, up to its matrix rank.
    Parameters
    ----------
    C : 2D array
        The covariance matrix.
    k : positive int, optional
        The number of eigenvectors to compute. ``None`` is interpreted as unspecified and
        yields all eigenvalues (up to the matrix rank of ``C``). The default is ``None``.
    method : non-negative int, optional
        Library subroutine to be utilized for eigenvector decomposition (EVD).
        The default is ``0``, unless ``C`` is ill-conditioned, then the default is ``1``.
        ``0``: Fastest method (incomplete EVD from ``scipy.linalg``)
        ``1``: Medium speed method (SVD from ``numpy.linalg``)
        ``2``: Medium speed method (SVD from ``scipy.linalg``)
        ``3``: Slowest method (complete EVD from ``numpy.linalg``)
        ``4``: Slowest method (complete EVD from ``scipy.linalg``)
    Returns
    -------
    U : 2D array
        The eigenvectors of the covariance matrix ``C``.
    S : 1D array
        The eigenvalues of the covariance matrix ``C``.
    '''

    # if k > rank(C)
    r = np.linalg.matrix_rank(C)
    if k == None or k > r:
        k = r

    if method == 0:
        if np.linalg.cond(C) > 1e7:
            # Default to method 1 if ill-conditioned
            method = 1
    
    ## Fastest method (incomplete EVD from scipy linalg)
    if method == 0:
        S, U = linalg.eigh(C, eigvals=(r - k, r - 1))
        # Sort eigenvalues and eigenvectors to non-increasing order
        S = S[::-1]
        U = U[:, ::-1]

    ## Medium speed methods (SVD)
    elif method == 1:
        U, S, _ = np.linalg.svd(C, full_matrices=False)

    elif method == 2:
        U, S, _ = linalg.svd(C, full_matrices=False)

    ## Slowest methods (complete EVD)
    elif method == 3:
        S, U = np.linalg.eig(C)
        # Sort eigenvalues and eigenvectors to non-increasing order
        ix = S.argsort()[::-1]
        S = S[ix]
        U = U[:, ix]

    elif method == 4:
        S, U = linalg.eig(C)
        # Take only real part of eigenvalues
        S = S.real

    if method != 0:
        U = U[:,:k]
        S = S[:k]

    return U, S


def normalize(data, axis=-2):
    ''' Normalizes columns/rows of data to unit L2 norm.
    '''
    norms = np.linalg.norm(data, axis=axis)

    if axis == -2:
        data = data * (1 / norms)
    elif axis == -1:
        data = data * (1 / norms[:, None])

    return data
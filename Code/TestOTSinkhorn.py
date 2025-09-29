# Author: Remi Flamary <remi.flamary@unice.fr>
#
# License: MIT License
# sphinx_gallery_thumbnail_number = 3

import numpy as np
import matplotlib.pylab as pl
import ot
import ot.plot
from ot.datasets import make_1D_gauss as gauss


def normalize_cost_matrix(M):
    M_min = np.min(M)
    M_max = np.max(M)
    M_normalized = (M - M_min) / (M_max - M_min)
    return M_normalized


def testSinkhorn(a, b, M, lambd):
    n = a.shape[0]
    m = b.shape[0]

    # Add a small value to avoid division by zero
    a = a.reshape(n, ) + 1e-10
    b = b.reshape(m, ) + 1e-10

    # Normalize the cost matrix
    M = normalize_cost_matrix(M)

    print(a.shape)
    print(b.shape)
    print(M.shape)

    Gs = ot.sinkhorn(a, b, M, lambd, verbose=True, numItermax=1000, stopThr=1e-9)
    print(Gs.shape)

    return Gs



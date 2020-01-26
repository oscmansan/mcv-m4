import numpy as np

import utils as h
import reconstruction as rc
import maths as mth
import fundamental as fd

def estimate_aff_hom(cams, vps):
    # your code here

    return aff_hom

def estimate_euc_hom(cams, vps):
    # your code here

    # make points homogeneous
    vpsh = [fd.make_homogeneous(vp) for vp in vps]

    # build A
    u = vpsh[0]
    v = vpsh[1]
    z = vpsh[2]
    A = [
    	[ u[0]*v[0], u[0]*v[1] + u[1]*v[0], u[0]*v[2] + u[2]*v[0], u[1]*v[1], u[1]*v[2] + u[2]*v[1], u[2]*v[2] ],
    	[ u[0]*z[0], u[0]*z[1] + u[1]*z[0], u[0]*z[2] + u[2]*z[0], u[1]*z[1], u[1]*z[2] + u[2]*z[1], u[2]*z[2] ],
    	[ v[0]*z[0], v[0]*z[1] + v[1]*z[0], v[0]*z[2] + v[2]*z[0], v[1]*z[1], v[1]*z[2] + v[2]*z[1], v[2]*z[2] ],
    	[ 0, 1, 0, 0, 0, 0 ],
    	[ 1, 0, 0, -1, 0, 0 ]

    ]

    # find w_v
    w_v = mth.nullspace(A)

    # build w
    w = [
    	[ w_v[0], w_v[1], w_v[2] ],
    	[ w_v[1], w_v[3], w_v[4] ],
    	[ w_v[2], w_v[4], w_v[5] ]
    ]

    # obtain K
    K = np.linalg.inv(np.linalg.cholesky(w))

    # obtain A
    M = cams[:,:3]
    A = np.linalg.cholesky(np.linalg.inv(M.T@w@M))

    # build euc_hom
    euc_hom = np.concatenate((np.concatenate((np.linalg.inv(A),np.zeros((1,3))), axis=0),np.asarray([[0, 0, 0, 1]]).T), axis=1)

    return euc_hom 

import cv2
import numpy as np

import utils as h
import maths as mth


def normalise_coord(p1, p2):
    # normalise both sets

    #alternative way of computing S1
    #dist_mean = np.mean(np.sqrt(((p1[:,:2] - mean_1)**2).sum(axis=1)), dtype=np.float32)
    #s1 = np.sqrt(2.)/dist_mean

    mean_1 = np.mean(p1[:,:2], axis=0, dtype=np.float32)
    S1 = np.sqrt(2.) / np.std(p1[:, :2], dtype=np.float32)
    T1 = np.float32(np.array([[S1, 0, -S1*mean_1[0]], [0, S1, -S1*mean_1[1]], [0, 0, 1]]))
    p1 = T1@p1.T

    mean_2 = np.mean(p2[:,:2],axis=0, dtype=np.float32)
    S2 = np.sqrt(2.) / np.std(p2[:, :2], dtype=np.float32)
    T2 = np.float32(np.array([[S2, 0, -S2*mean_2[0]], [0, S2, -S2*mean_2[1]], [0, 0, 1]]))
    p2 = T2@p2.T

    if h.debug >= 0:
        print("    Coordinates normalised")

    return p1.T, p2.T, T1, T2


def compute_fundamental(p1, p2, eight_alg):
    # compute fundamental matrix with normalised coordinates

    tol_rsc = np.array([1.5, 1.5, 1])

    if h.normalise:
        # make coordinates homogeneous:
        p1h = make_homogeneous(p1)
        p2h = make_homogeneous(p2)
        # normalise coordinates 
        p1_norm, p2_norm, T1, T2 = normalise_coord(p1h, p2h)
        #tolerance normalised for RANSAC method (not used by LMEDS, left for
        #compatibility with normalise=False path)
        tol_rsc_nrm = T1@tol_rsc
        # Only LMEDS seems to work well with normalised coordinates
        fund_method = cv2.FM_LMEDS
    else:
        p1_norm, p2_norm = p1, p2
        tol_rsc_nrm = tol_rsc
        fund_method = cv2.FM_RANSAC

    if eight_alg:
        fund_method = cv2.FM_8POINT

    F, mask = cv2.findFundamentalMat(p1_norm[:, :2], p2_norm[:, :2], fund_method, tol_rsc_nrm[0], 0.99)

    if h.normalise:
        # denormalise F
        F = T2.T@F@T1 
        F = F / F[2][2]
  
    if h.debug >= 0:
        print('    Fundamental Matrix estimated')
    if h.debug > 1:
        print("      Fundamental Matrix: \n", F)

    return F, mask


def apply_mask(x1, x2, mask, F):
    # use F mask for filtering out outliers
    if h.debug > 2: 
        print("before F mask:\n")
        xh1 = make_homogeneous(x1)
        xh2 = make_homogeneous(x2)
        mth.print_epipolar_eq(xh1, xh2, F)

    # apply mask of inliers to the set of matches
    x1 = x1[mask.ravel() == 1]
    x2 = x2[mask.ravel() == 1]

    if h.debug > 2: 
        print("after F mask:\n")
        xh1 = make_homogeneous(x1)
        xh2 = make_homogeneous(x2)
        mth.print_epipolar_eq(xh1, xh2, F)

    if h.debug >= 0:
        print("    Mask given by F applied")
    if h.debug > 0:
        print("      F mask has selected", x1.shape[0], "inliers")

    return x1, x2


def refine_matches(x1, x2, F):
    # use the optimal triangulation method (Algorithm 12.1 from MVG)
    nx1, nx2 = cv2.correctMatches(F, np.reshape(x1, (1, -1, 2)), np.reshape(x2, (1, -1, 2)))

    # get the points back in matrix configuration
    xr1 = np.float32(np.reshape(nx1,(-1, 2)))
    xr2 = np.float32(np.reshape(nx2,(-1, 2)))

    if h.debug >= 0:
        print("  Matches corrected with Optimal Triangulation Method")

    if h.debug > 2: 
        print("xr1: \n", xr1)
        print("xr2: \n", xr2)
        print("after correctMatches: ")
        xrh1 = make_homogeneous(xr1)
        xrh2 = make_homogeneous(xr2)
        mth.print_epipolar_eq(xrh1, xrh2, F)

    return xr1.T, xr2.T 


def search_more_matches(out1, out2, F):
    # your code here

    e = 0.00155

    outh1 = make_homogeneous(out1)
    outh2 = make_homogeneous(out2)

    xn1 = np.empty([0, 2], dtype=np.int32)
    xn2 = np.empty([0, 2], dtype=np.int32)
    on1 = np.empty([0, 2], dtype=np.int32)
    on2 = np.empty([0, 2], dtype=np.int32)

    for oh1, oh2 in zip(outh1, outh2):
        # compute epipolar lines
        l1 = F.T@oh2
        l2 = F@oh1

        # distance from a point to a line
        d1 = abs(np.dot(l1, oh1)) / np.sqrt(np.sum(l1[0:2]**2))
        d2 = abs(np.dot(l2, oh2)) / np.sqrt(np.sum(l2[0:2]**2))
        d = d1 + d2

        if d < e:
            xn1 = np.r_[xn1, [np.int32(oh1[0:2]/oh1[2])]]
            xn2 = np.r_[xn2, [np.int32(oh2[0:2]/oh2[2])]]
        else:
            on1 = np.r_[on1, [np.int32(oh1[0:2]/oh1[2])]]
            on2 = np.r_[on2, [np.int32(oh2[0:2]/oh2[2])]]

    return xn1, xn2, on1, on2


def make_homogeneous(p):
    if p.shape[1] != 2:
        print("WARNING - Coordinates have ", p.shape, " dimensions: not made homogeneous")
        return p
    
    # Add homogeneous coordinates 
    hom_c = np.ones_like(p[:, 1], dtype=np.float32)
    p = np.c_[p, hom_c]

    # Perhaps is quicker to use opencv implementation, but it's 
    # cumbersome... seems to use lists
    #p = np.reshape(cv2.convertPointsToHomogeneous(p), (-1, 3))

    if h.debug > 0:
        print("    Coordinates made homogeneous")

    return p

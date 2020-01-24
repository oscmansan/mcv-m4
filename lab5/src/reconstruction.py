import cv2
import numpy as np

import utils as h
import maths as mth

def compute_proj_camera(F, i):
    # Result 9.15 of MVG (v = 0, lambda = 1). It assumes P1 = [I|0]
    # your code here

    # compute epipole
    e = mth.nullspace(F.T)

    # build [e]_x
    e_x = mth.hat_operator(e)

    # compute P
    P = e_x@F
    P = np.concatenate((P, e), axis=1)

    return P

def estimate_3d_points(P1, P2, xr1, xr2):
    # Triangulate 3D points from camera matrices
    Xprj = cv2.triangulatePoints(P1, P2, xr1, xr2) 

    # Divide by the last column 
    Xprj = Xprj / Xprj[3, :]

    if h.debug >2:
        print("  X estimated:\n", Xprj)

    return Xprj

def compute_reproj_error(X, x1, x2, cam): 
    # your code here

    # initialize variables
    error = 0

    # reshape X
    X = np.reshape(X, (X.shape[1], X.shape[0]))

    for i, X_point in enumerate(X):
        # compute reprojected points
        x1_reprojected = cam[0]@X_point
        x2_reprojected = cam[1]@X_point
        # convert to euclidean coordinates
        x1_reprojected = [x1_reprojected[0] / x1_reprojected[2], x1_reprojected[1] / x1_reprojected[2]]
        x2_reprojected = [x2_reprojected[0] / x2_reprojected[2], x2_reprojected[1] / x2_reprojected[2]]
        # compute reprojection error
        e1 = (x1[i][0] - x1_reprojected[0]) ** 2 + (x1[i][1] - x1_reprojected[1]) ** 2
        e2 = (x2[i][0] - x2_reprojected[0]) ** 2 + (x2[i][1] - x2_reprojected[1]) ** 2
        e = e1 + e2
        # accumulate
        error += e

    return error

def transform(aff_hom, Xprj, cams_pr):
    # your code here

    return Xaff, cams_aff

def resection(tracks, img):
    # your code here

    return P

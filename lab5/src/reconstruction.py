import cv2
import numpy as np

import utils as h
import maths as mth

def compute_proj_camera(F, i):
    # Result 9.15 of MVG (v = 0, lambda = 1). It assumes P1 = [I|0]
    # your code here

    # compute epipole e'
    e = mth.nullspace(F.T)

    # build [e']_x
    ske = mth.hat_operator(e)

    # compute P
    P = np.concatenate((ske@F, e), axis=1)

    return P

def estimate_3d_points(P1, P2, xr1, xr2):
    # Triangulate 3D points from camera matrices
    Xprj = cv2.triangulatePoints(P1, P2, xr1, xr2) 

    # Divide by the last column 
    Xprj = Xprj / Xprj[3, :]

    if h.debug >2:
        print("  X estimated:\n", Xprj)

    return Xprj

def compute_reproj_error(X, P1, P2, xr1, xr2):
    # your code here

    xp1 = P1@X
    xp2 = P2@X
    xp1 = xp1[0:2, :] / xp1[None, 2, :]
    xp2 = xp2[0:2, :] / xp2[None, 2, :]

    error = np.sum(np.sum((xr1-xp1)**2)+np.sum((xr2-xp2)**2))

    return error

def transform(aff_hom, Xprj, cams_pr):
    # your code here

    return Xaff, cams_aff

def resection(tracks, img):
    # your code here

    return P

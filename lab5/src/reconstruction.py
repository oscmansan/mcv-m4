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
    e_x = np.asarray([[0, -e[2], e[1]],
           [e[2], 0, -e[0]],
           [-e[1], e[0], 0]
    ])

    # compute P
    P = e_x@F
    P = np.concatenate((P,e),axis=1)

    return P

def estimate_3d_points(P1, P2, xr1, xr2):
    # Triangulate 3D points from camera matrices
    Xprj = cv2.triangulatePoints(P1, P2, xr1, xr2) 

    # Divide by the last column 
    Xprj = Xprj / Xprj[3, :]

    if h.debug >2:
        print("  X estimated:\n", Xprj)

    return Xprj

def compute_reproj_error(X, cam): 
    # your code here 

    return error

def transform(aff_hom, Xprj, cams_pr):
    # your code here

    return Xaff, cams_aff

def resection(tracks, img):
    # your code here

    return P

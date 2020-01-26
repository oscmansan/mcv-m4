import cv2
import numpy as np

import utils as h
import maths as mth


def compute_proj_camera(F, i):
    # Result 9.15 of MVG (v = 0, lambda = 1). It assumes P1 = [I|0]

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
    # project 3D points using P
    xp1 = P1@X
    xp2 = P2@X
    xp1 = euclid(xp1.T).T
    xp2 = euclid(xp2.T).T

    # compute reprojection error
    error = np.sum(np.sum((xr1-xp1)**2)+np.sum((xr2-xp2)**2))

    return error


def transform(aff_hom, Xprj, cams_pr):
    # Algorithm 19.2 of MVG

    Xaff = aff_hom@Xprj
    Xaff = Xaff / Xaff[3, :]

    cams_aff = [cam@np.linalg.inv(aff_hom) for cam in cams_pr]

    return Xaff, cams_aff


def resection(tracks, i):
    # extract 3D-2D correspondences from tracks
    pts3d = []
    pts2d = []
    for tk in tracks:
        if i in tk.views:
            pts3d.append(tk.pt)
            pts2d.append(tk.views[i])
    pts3d = np.array(pts3d)
    pts2d = np.array(pts2d)

    # convert to homogeneous coordinates
    pts3d = homog(pts3d)
    pts2d = homog(pts2d)

    # normalize points
    pts3d, T1 = normalize3dpts(pts3d)
    pts2d, T2 = normalize2dpts(pts2d)

    # DLT algorithm
    n = 6  # minimal solution
    A = np.empty((2*n, 12))
    for i in range(n):
        X = pts3d[i]
        x, y, w = pts2d[i]
        A[2*i, :] = np.concatenate((np.zeros(4), -w*X, y*X))
        A[2*i+1, :] = np.concatenate((w*X, np.zeros(4), -x*X))

    U, D, V = np.linalg.svd(A)
    p = V[:, -1]
    P = p.reshape((3, 4))

    # TODO: over-determine solution (n>6)
    # TODO: minimize geometric error
    # Note: RANSAC is not needed since points in tracks are inliers

    # denormalize P
    P = T2.T @ P @ T1
    P /= P[-1, -1]

    if h.debug >= 0:
        print('    Camera Matrix estimated')
    if h.debug > 1:
        print('      Camera Matrix: {}\n'.format(P))

    return P


def homog(x):
    return np.concatenate((x, np.ones((x.shape[0], 1))), axis=1)


def euclid(x):
    return x[:, :-1] / x[:, [-1]]


def normalize2dpts(pts):
    mean = np.mean(pts[:, :2], axis=0, dtype=np.float32)
    S = np.sqrt(2.) / np.std(pts[:, :2], dtype=np.float32)
    T = np.float32(np.array([[S, 0, -S * mean[0]],
                             [0, S, -S * mean[1]],
                             [0, 0, 1]]))
    pts = T @ pts.T
    return pts.T, T


def normalize3dpts(pts):
    mean = np.mean(pts[:, :3], axis=0, dtype=np.float32)
    S = np.sqrt(2.) / np.std(pts[:, :3], dtype=np.float32)
    T = np.float32(np.array([[S, 0, 0, -S * mean[0]],
                             [0, S, 0, -S * mean[1]],
                             [0, 0, S, -S * mean[2]],
                             [0, 0, 0, 1]]))
    pts = T @ pts.T
    return pts.T, T

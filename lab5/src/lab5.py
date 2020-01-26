# Week 5. 3D reconstruction from N non calibrated cameras. 
# This exercise implements a vanilla version of Structure from Motion technique. 

# Use python3 to execute this programme: 
#    $ python3 lab5.py <number of images to process>

import sys

# opencv, numpy
import cv2
import numpy as np

# bundle adjustment
import pysba as ba

# project files
import utils as h
import maths as mth
import matches as mt
import fundamental as fd
import track as tk
import vps as vp
import autocalibration as ac
import reconstruction as rc

def main(argv):

    # usage
    if (len(argv) != 2):
        print ("Usage: python3 lab5.py <number of images to process>")
        sys.exit(0)

    imgs = []        # list of images
    feats = []       # list of features
    matches = []     # list of dictionaries
    tracks = []      # list of tracking views 
    hs_vs = {}       # dictionary as hash table of views-tracks
    vps = []         # list of vanishing points 
    cams_pr = []     # list of projective cameras
    cams_aff = []    # list of affine cameras
    cams_euc = []    # list of euclidean cameras
    Xpr = []         # list of projective 3d points
    Xaff = []        # list of affine 3d points
    Xeuc = []        # list of euclidean 3d points

    # Get number of images to process
    n = int(argv[1])
    for i in range(0,n):
        if h.debug >=0:
            print("Processing image", i, "of sequence")

        # read image
        imgs.append(h.read_image(i))

        # find features
        feati = mt.find_features(imgs[i], i)
        feats.append(feati)

        if i == 0:
            P0 = np.float32(np.c_[np.eye(3), np.zeros(3)])
            cams_pr.append(P0)
            vps.append(vp.estimate_vps(imgs[i]))
            if h.debug >= 0:
                print ("  Camera 0 set to identity")
        else:
            for prev in range(len(imgs)-2, -1, -1): # from the penultimate image
                if h.debug >=0:
                    print("  Matching images",prev, "and",i,"for obtaining tracks")
                # match features
                m_ij = mt.match_features(feats[prev][1], feats[i][1], prev, i)
                m_ijf = mt.filter_matches(feats[prev][0], feats[i][0], m_ij, prev, i)

                # Iterative process for refining F and detecting all good matches
                if h.debug >= 0:
                    print ("  Using epipolar constraint to find more matches")
                # Prepare variables to be used in the iteration
                incr_match = 1   
                eight_alg = False
                # inliers
                x1 = m_ijf[0]
                x2 = m_ijf[1]
                # outliers
                o1 = m_ijf[2]
                o2 = m_ijf[3]

                while incr_match > 0:
                    # find fundamental matrix
                    F, mask = fd.compute_fundamental(x1, x2, eight_alg)
                    # apply mask
                    xs1, xs2 = fd.apply_mask(x1, x2, mask, F)
                    # look for new matches
                    # TODO implement a search for new inliers based on epipolar
                    # error given by F (MVG, Alg 11.4 (v)). Set the maximum error that an inlier may
                    # throw on the epipolar equation to 0.00155
                    # Returns: 
                    #   xn1, xn2: new inliers
                    #   o1, o2: still outliers
                    xn1, xn2, o1, o2 = fd.search_more_matches(o1, o2, F)

                    # join the new matches to the inliers
                    x1 = np.concatenate((xs1, xn1), axis=0)
                    x2 = np.concatenate((xs2, xn2), axis=0)
                    incr_match = xn1.shape[0]
                    eight_alg = True

                if h.debug >= 0:
                    print ("  Search with epipolar constraint finished")
                if h.debug >0:
                    print("    Inliers have grown to", x1.shape[0])

                # refine matches and update the contents of matches
                xr1, xr2 = fd.refine_matches(x1, x2, F)
                tk.add_tracks(x1, x2, xr1.T, xr2.T, prev, i, tracks, hs_vs)

                if h.debug >=0:
                    print("  Tracks added after matching",prev,"and",i)
                if h.debug >0:
                    print("    Size of tracks:",len(tracks))
                    print("    Size of hash table of views:",len(hs_vs))

                if h.debug_display:
                    h.display_epilines(imgs[prev], imgs[i], x1, x2, F)
                    h.show_matches(imgs[prev], imgs[i], x1, x2)

            # compute projective cameras to use in projective reconstruction
            if i == 1:
                # TODO Compute the projective camera given F, according to
                # Result 9.15 of MVG (v = 0, lambda = 1).
                cams_pr.append(rc.compute_proj_camera(F, i))
            else:
                # TODO Compute resection as in MVG, Alg 7.1
                cams_pr.append(rc.resection(tracks, i))
            if h.debug >= 0:
                print("  Resection of camera", i, "performed")

            # projective triangulation for 3D structure
            Xprj = rc.estimate_3d_points(cams_pr[i-1], cams_pr[i], xr1, xr2)
            if h.debug >= 0:
                print('  Projective reconstruction estimated')

            # TODO Add estimated 3d projective points to tracks
            tk.add_pts_tracks(Xprj, x1, x2, tracks, hs_vs)
            if h.debug >= 0:
                print('  Projective 3D points added to tracks')

            # TODO compute projective reprojection error
            error_prj = rc.compute_reproj_error(Xprj, cams_pr[i-1], cams_pr[i], xr1, xr2)
            if h.debug >0:
                print("    Projective reprojection error:", error_prj)
            
            # Affine rectification
            vps.append(vp.estimate_vps(imgs[i]))
            # TODO Estimate homografy that makes an affine rectification
            # With the vanishing points, the plane at the infinity is computed. 
            # Then the affine homography is built with the coordinates of the infinity plane
            aff_hom = ac.estimate_aff_hom(cams_pr[i-1:], vps[i-1:])

            # TODO Transform 3D points and cameras to affine space
            Xaff, cams_aff = rc.transform(aff_hom, Xprj, cams_pr)

            # TODO Add estimated 3d affine points to tracks (reuse your code)
            tk.add_pts_tracks(Xaff, x1, x2, tracks, hs_vs)
            if h.debug >= 0:
                print('  Affine 3D points added to tracks')
            
            # TODO compute affine reprojection error (reuse your code)
            error_aff = rc.compute_reproj_error(Xaff, cams_aff[i-1], cams_aff[i], xr1, xr2)
            if h.debug >0:
                print("    Affine reprojection error:", error_aff)

            # Metric rectification
            # TODO Perform Metric rectification. First compute the transforming
            # homography from vanishing points and the camera constrains skew = 0,
            # squared pixels. Then perform the transformation to Euclidean space
            # (reuse your code)
            euc_hom = ac.estimate_euc_hom(cams_aff[i], vps[i])
            Xeuc, cams_euc = rc.transform(euc_hom, Xaff, cams_aff)

            # TODO Add estimated 3d euclidean points to tracks (reuse your code)
            tk.add_pts_tracks(Xeuc, x1, x2, tracks, hs_vs)
            if h.debug >= 0:
                print('  Euclidean 3D points added to tracks')
            
            # TODO compute metric reprojection error (reuse your code)
            error_euc = rc.compute_reproj_error(Xeuc, cams_euc[i-1], cams_euc[i], xr1, xr2)
            if h.debug >0:
                print("    Euclidean reprojection error:", error_euc)

            # Bundle Adjustment
            # TODO Adapt cameras and 3D points to PySBA format
            cams_ba, X_ba, x_ba, cam_idxs, x_idxs = ba.adapt_format_pysba(tracks)
            badj = ba.PySBA(cams_ba, X_ba, x_ba, cam_idxs, x_idxs)
            cams_euc, Xeuc = badj.bundleAdjust()
            # TODO Update 3D points and tracks with optimised cameras and points
            tk.update_ba_pts_tracks(Xeuc, tracks)
            if h.debug >=0:
                print("  Bundle Adjustment performed over", i,"images")

            # render results
            if h.debug_display:
                h.display_3d_points(Xeuc.T[:,:3])

    if h.debug >=0:
        print ("Structure from Motion applied on sequence of", n, "images")

"""
Optional tasks
    - Estimate affine homography from the 3 vanishing points and F (Alg. 13.1
      p332, result 10.3 p271)
    - Perform Bundle Adjustment over the estimation of the vanishing points and
      all available images, with PySBA
    - Implement a more sophisticated resection method: 
        -Xiao-Shan Gao, Xiao-Rong Hou, Jianliang Tang, and Hang-Fei Cheng. 
         Complete solution classification for the perspective-three-point problem. 
         Pattern Analysis and Machine Intelligence, IEEE Transactions on, 25(8):930–943, 2003
        -Tong Ke and Stergios Roumeliotis. An efficient algebraic solution to the perspective-three-point 
         problem. In Computer Vision and Pattern Recognition (CVPR), 2017 IEEE Conference on. IEEE, 2017
    - Investigate better implementations for the "tracks" structure
    - Investigate strategies to improve the pipeline:
            -in results: number of points, reprojection error, camera poses
            -in implementation: time of computation, resources, etc.
      Reference papers for improvement strategies: 
            -J. L. Schönberger and J. Frahm, "Structure-from-Motion Revisited," 2016 IEEE Conference 
             on Computer Vision and Pattern Recognition (CVPR), Las Vegas, NV, 2016, pp. 4104-4113.
"""

main(sys.argv)

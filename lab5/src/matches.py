import cv2
import numpy as np

import utils as h

def find_features(img, i):
    # find the keypoints and descriptors in img
    #SIFT used
    sift = cv2.xfeatures2d.SIFT_create()

    kp, des = sift.detectAndCompute(img, None)

    if h.debug >= 0:
        print ("  Features detected in image",i)
    if (h.debug > 0):
        print("    Found", len(kp), "features ")

    return kp, des # this is a tuple

def match_features(des1, des2, i, j):
    # FLANN parameters
    FLANN_INDEX_KDTREE = 0
    index_params = dict(algorithm = FLANN_INDEX_KDTREE, trees = 5)
    search_params = dict(checks=50)

    flann = cv2.FlannBasedMatcher(index_params,search_params)
    matches = flann.knnMatch(des1, des2, k=2)

    if h.debug >= 0:
        print ("  Correspondences matched between images", i, "and", j)
    if (h.debug > 0):
        print ("    Found", len(matches), "matching correspondences")

    return matches

def filter_matches(kp1, kp2, matches, imgi, imgj):
    x1 = np.empty([0, 2], dtype=np.int32)
    x2 = np.empty([0, 2], dtype=np.int32)
    o1 = np.empty([0, 2], dtype=np.int32)
    o2 = np.empty([0, 2], dtype=np.int32)

    # ratio test as per Lowe's paper
    for m, n in matches: 
        if m.distance < 0.8*n.distance:
            x1 = np.r_[x1, [np.int32(np.array(kp1[m.queryIdx].pt))]]
            x2 = np.r_[x2, [np.int32(np.array(kp2[m.trainIdx].pt))]]
        else:
            o1 = np.r_[o1, [np.int32(np.array(kp1[m.queryIdx].pt))]]
            o2 = np.r_[o2, [np.int32(np.array(kp2[m.trainIdx].pt))]]

    if h.debug >= 0:
        print ("  Matches between", imgi, "and", imgj, "filtered with Lowe's ratio")
    if (h.debug > 0):
        print ("    Selected", x1.shape[0], "matches")

    return [x1, x2, o1, o2]


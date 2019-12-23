function [mosaic] = compute_mosaic_3(ima_path, imb_path, imc_path, threshold, corners)

% Read images
imargb = imread(ima_path);
imbrgb = imread(imb_path);
imcrgb = imread(imc_path);

% Convert to grayscale if necessary
if size(imargb,3) == 3
    ima = sum(double(imargb), 3) / 3 / 255;
else
    ima = double(imargb);
end
if size(imbrgb,3) == 3
    imb = sum(double(imbrgb), 3) / 3 / 255;
else
    imb = double(imbrgb);
end
if size(imcrgb,3) == 3
    imc = sum(double(imcrgb), 3) / 3 / 255;
else
    imc = double(imcrgb);
end
    
% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

% Match SIFT keypoints
matches_ab = siftmatch(desc_a, desc_b);
matches_bc = siftmatch(desc_b, desc_c);

% Compute homography
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, ~] = ransac_homography_adaptive_loop(xab_a, xab_b, threshold, 1000);
xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, ~] = ransac_homography_adaptive_loop(xbc_b, xbc_c, threshold, 1000);

% Build mosaic
iwb = apply_H_v2(imbrgb, eye(3), corners);
iwa = apply_H_v2(imargb, Hab, corners);
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);

mosaic = max(iwc, max(iwb, iwa));

end
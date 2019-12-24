function [mosaic] = compute_mosaic_3(images, threshold, corners)

% Read images
imargb = imread(images{1});
imbrgb = imread(images{2});
imcrgb = imread(images{3});

% Convert to grayscale if necessary
if size(imargb,3) == 3
    ima = rgb2gray(imargb);
else
    ima = imargb;
end
if size(imbrgb,3) == 3
    imb = rgb2gray(imbrgb);
else
    imb = imbrgb;
end
if size(imcrgb,3) == 3
    imc = rgb2gray(imcrgb);
else
    imc = imcrgb;
end

ima = double(ima) / 255;
imb = double(imb) / 255;
imc = double(imc) / 255;
    
% Compute SIFT keypoints
fprintf('Computing SIFT keypoints...\n');
tic;
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);
toc;

% Match SIFT keypoints
fprintf('Matching SIFT keypoints...\n');
tic;
matches_ab = siftmatch(desc_a, desc_b);
matches_bc = siftmatch(desc_b, desc_c);
toc;

% Compute homographies
fprintf('Computing homographies...\n');
tic;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, ~] = ransac_homography_adaptive_loop(xab_a, xab_b, threshold, 1000);
xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, ~] = ransac_homography_adaptive_loop(xbc_b, xbc_c, threshold, 1000);
toc;

% Build mosaic
fprintf('Building mosaic...\n');
tic;
iwb = apply_H_v2(imbrgb, eye(3), corners);
iwa = apply_H_v2(imargb, Hab, corners);
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);
toc;

mosaic = max(iwc, max(iwb, iwa));

end
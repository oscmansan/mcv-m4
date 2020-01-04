function H = dlt(ima,imb)
%DLT Summary of this function goes here
%   Detailed explanation goes here

% convert to grayscale if necessary
if size(ima,3) == 3
    ima = rgb2gray(ima);
end
if size(imb,3) == 3
    imb = rgb2gray(imb);
end

ima = double(ima) / 255;
imb = double(imb) / 255;

% compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01, 'Verbosity', 1);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01, 'Verbosity', 1);

% match SIFT keypoints
matches = siftmatch(desc_a, desc_b);

% compute homography (normalized DLT)
th = 3;
x_a = [points_a(1:2, matches(1,:)); ones(1, length(matches))];
x_b = [points_b(1:2, matches(2,:)); ones(1, length(matches))];
[H, ~] = ransac_homography_adaptive_loop(x_a, x_b, th, 1000);

end

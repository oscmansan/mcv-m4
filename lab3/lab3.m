%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 3: The geometry of two views 
% (application: photo-sequencing)

addpath('../lab2/sift'); % ToDo: change 'sift' to the correct path where you have the sift functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute the fundamental matrix

% Two camera matrices for testing purposes
P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = P1 * X;
x2_test = P2 * X;

% Estimated fundamental matrix
% ToDo: create the following function that estimates F using the normalised 8 point algorithm
F_es = fundamental_matrix(x1_test, x2_test);

% Real fundamental matrix
% ToDo: write the expression of the real fundamental matrix for P1 and P2
F_gt = [0, -t(3), t(2);
        t(3), 0, -t(1);
        -t(2), t(1), 0] * R;  % inv(Kp).' * [tp]_x * R * inv(K), where K = eye(3,3)

% Evaluation: these two matrices should be very similar
F_es = F_es / norm(F_es);
F_gt = F_gt / norm(F_gt);

% Since F is up-to-scale, the sign might be inverted
assert(iscloseall(F_es,F_gt) || iscloseall(-F_es,F_gt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Robustly fit fundamental matrix

% Read images
im1rgb = imread('Data/0000_s.png');
im2rgb = imread('Data/0001_s.png');
im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;

% show images
figure;
subplot(1,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(1,2,2); imshow(im2rgb); axis image; title('Image 2');


%% Compute SIFT keypoints
% (if you have problems with SIFT you may use SURF features.
% with Matlab functions 'detectSURFFeatures', 'extractFeatures' and 'matchFeatures'
% more details in https://ch.mathworks.com/help/vision/ref/matchfeatures.html)

% (make sure that the sift folder provided in lab2 is on the path)

[points_1, desc_1] = sift(im1, 'Threshold', 0.01);
[points_2, desc_2] = sift(im2, 'Threshold', 0.01);

%% Match SIFT keypoints between a and b
matches = siftmatch(desc_1, desc_2);
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches, 'Stacking', 'v');

% p1 and p2 contain the homogeneous coordinates of the matches
p1 = [points_1(1:2, matches(1,:)); ones(1, length(matches))];
p2 = [points_2(1:2, matches(2,:)); ones(1, length(matches))];

% ToDo: create this function (you can use as a basis 'ransac_homography_adaptive_loop.m')
[F, inliers] = ransac_fundamental_matrix(p1, p2, 2.0);

% show inliers
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches(:,inliers), 'Stacking', 'v');
title('Inliers');

vgg_gui_F(im1rgb, im2rgb, F');


%% Plot some epipolar lines

l2 = F*p1; % epipolar lines in image 2 % ToDo
l1 = F.'*p2; % epipolar lines in image 1 % ToDo

% choose three random indices
m1 = inliers(10);
m2 = inliers(20);
m3 = inliers(30);

% image 1 (plot the three points and their corresponding epipolar lines)
figure;
imshow(im1rgb);
hold on;
plot(p1(1, m1), p1(2, m1), '+g');
plot_homog_line(l1(:, m1));

plot(p1(1, m2), p1(2, m2), '+g');
plot_homog_line(l1(:, m2));

plot(p1(1, m3), p1(2, m3), '+g');
plot_homog_line(l1(:, m3));

% image 2 (plot the three points and their corresponding epipolar lines)
figure;
imshow(im2rgb);
hold on;
plot(p2(1, m1), p2(2, m1), '+g');
plot_homog_line(l2(:, m1));

plot(p2(1, m2), p2(2, m2), '+g');
plot_homog_line(l2(:, m2));

plot(p2(1, m3), p2(2, m3), '+g');
plot_homog_line(l2(:, m3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Photo-sequencing with aerial images

% In this part we will compute a simplified version of the algorithm
% explained in the Photo-sequencing paper. 
% Since we do not have two images
% taken from roughly the same viewpoint at two different time instants we
% will manually pick a dynamic point corresponding to a point in a van 
% (identified by index 'idx_car_I1') and the projection of its 3D trajectory 
% in the reference image. Then we will compute the projection (to the reference image) 
% of three points on this 3D trajectory at three different time instants 
% (corresponding to the time when the three other provided images where taken). 

clear all;

% Read images
im1rgb = imread('Data/frame_00000.tif');
im2rgb = imread('Data/frame_00001.tif');
im3rgb = imread('Data/frame_00002.tif');
im4rgb = imread('Data/frame_00003.tif');

im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;

% show images
figure;
subplot(2,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(2,2,2); imshow(im2rgb); axis image; title('Image 2');
subplot(2,2,3); imshow(im3rgb); axis image; title('Image 3');
subplot(2,2,4); imshow(im4rgb); axis image; title('Image 4');

% Compute SIFT keypoints
[points_1, desc_1] = sift(im1, 'Threshold', 0.015); % Do not change this threshold!
[points_2, desc_2] = sift(im2, 'Threshold', 0.015);
[points_3, desc_3] = sift(im3, 'Threshold', 0.015);
[points_4, desc_4] = sift(im4, 'Threshold', 0.015);

%% Take image im1 as reference image (image 1) and compute the fundamental 
% matrices needed for computing the trajectory of point idx_car_I1
% (use the SIFT keypoints previously computed)

% Compute sift matches
matches_12 = siftmatch(desc_1, desc_2);
matches_13 = siftmatch(desc_1, desc_3);
matches_14 = siftmatch(desc_1, desc_4);

% F matrix of image 2
p1 = [points_1(1:2, matches_12(1,:)); ones(1, length(matches_12))];
p2 = [points_2(1:2, matches_12(2,:)); ones(1, length(matches_12))];
% F2 = fundamental_matrix(p1, p2);
[F2, inliers_2] = ransac_fundamental_matrix(p1, p2, 2.0);

% F matrix of image 3
p1 = [points_1(1:2, matches_13(1,:)); ones(1, length(matches_13))];
p2 = [points_3(1:2, matches_13(2,:)); ones(1, length(matches_13))];
% F3 = fundamental_matrix(p1, p2);
[F3, inliers_3] = ransac_fundamental_matrix(p1, p2, 2.0);

% F matrix of image 4
p1 = [points_1(1:2, matches_14(1,:)); ones(1, length(matches_14))];
p2 = [points_4(1:2, matches_14(2,:)); ones(1, length(matches_14))];
% F4 = fundamental_matrix(p1, p2);
[F4, inliers_4] = ransac_fundamental_matrix(p1, p2, 2.0);

vgg_gui_F(im1rgb, im2rgb, F2');
vgg_gui_F(im1rgb, im3rgb, F3');
vgg_gui_F(im1rgb, im4rgb, F4');

%% Plot the car trajectory (keypoint idx_car_I1 in image 1)

% ToDo: complete the code

idx_car_I1 = 1197;
idx_car_I2 = matches_12(2, (find(matches_12(1,:) == idx_car_I1)));
idx_car_I3 = matches_13(2, (find(matches_13(1,:) == idx_car_I1)));
idx_car_I4 = matches_14(2, (find(matches_14(1,:) == idx_car_I1)));

% coordinates (in image 1) of the keypoint idx_car_I1 (point in a van). 
% point1_1 is the projection of a 3D point in the 3D trajectory of the van
point1_1 = [points_1(1:2,idx_car_I1)' 1]';
% coordinates (in image 1) of another 3D point in the same 3D trajectory of
% the van
point1_2 = [334 697 1]'; % (this is a given data)

% l1 is the projection of the 3D trajectory of keypoint idx_car_I1
% (it is the line that joins point1_1 and point1_2)
l1 = cross(point1_1, point1_2); % ToDo: compute the line
% plot the line
figure;imshow(im1);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(points_1(1,1197), points_1(2,1197), 'y*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 2
point2 = [points_2(1:2,idx_car_I2)' 1]';
% ToDo: compute the epipolar line of point2 in the reference image
l2 = F2.'*point2;
% plot the epipolar line
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'c');
% ToDo: compute the projection of point idx_car_I2 in the reference image 
pi2 = cross(l1,l2);
% plot this point
plot(pi2(1)/pi2(3), pi2(2)/pi2(3), 'c*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 3
point3 = [points_3(1:2,idx_car_I3)' 1]';
% ToDo: compute the epipolar line of point3 in the reference image
l3 = F3.'*point3;
% plot the epipolar line
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
% ToDo: compute the projection of point idx_car_I3 in the reference image
pi3 = cross(l1,l3);
plot(pi3(1)/pi3(3), pi3(2)/pi3(3), 'b*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 4
point4 = [points_4(1:2,idx_car_I4)' 1]';
% ToDo: compute the epipolar line of point4 in the reference image
l4 = F4.'*point4;
% plot the epipolar line
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
% ToDo: compute the projection of point idx_car_I4 in the reference image
pi4 = cross(l1,l4);
plot(pi4(1)/pi4(3), pi4(2)/pi4(3), 'g*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. OPTIONAL: Photo-sequencing with your own images

% 4.1 Take a set of images of a moving scene from different viewpoints at 
%     different time instants. At least two images have to be taken from
%     roughly the same location by the same camera.
%
% 4.2 Implement the first part (until line 16) of the Algorithm 1 of the 
%     Photo-sequencing paper with a selection of the detected dynamic
%     features. You may reuse the code generated for the previous question.
%

% Read images
im1rgb=imread('Data/SkateBoard/IMG0_006.png');
im2rgb=imread('Data/SkateBoard/IMG0_004.png');
im3rgb=imread('Data/SkateBoard/IMG0_003.png');
im4rgb=imread('Data/SkateBoard/IMG0_002.png');

% Convert to grayscale
im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;

figure;
subplot(2,2,1); imshow(im1); axis image; title('1st image');
subplot(2,2,2); imshow(im2); axis image; title('2nd image');
subplot(2,2,3); imshow(im3); axis image; title('3rd image');
subplot(2,2,4); imshow(im4); axis image; title('4rd image');

% Compute SIFT descriptors
[points_1, desc_1] = sift(im1, 'Threshold', 0.015);
[points_2, desc_2] = sift(im2, 'Threshold', 0.015);
[points_3, desc_3] = sift(im3, 'Threshold', 0.015);
[points_4, desc_4] = sift(im4, 'Threshold', 0.015);

% Match images
matches_12 = siftmatch(desc_1, desc_2);
matches_13 = siftmatch(desc_1, desc_3);
matches_14 = siftmatch(desc_1, desc_4);

%% Compute Euclidean distance between matched points
d_12 = vecnorm(points_1(1:2,matches_12(1,:))-points_2(1:2,matches_12(2,:)),2,1);
d_13 = vecnorm(points_1(1:2,matches_13(1,:))-points_3(1:2,matches_13(2,:)),2,1);
d_14 = vecnorm(points_1(1:2,matches_14(1,:))-points_4(1:2,matches_14(2,:)),2,1);

% Select static points
th = 20;  % chosen empirically (in pixels)
idx_static_1 = matches_12(1,d_12<th);
matches_static_12 = matches_12(:,ismember(matches_12(1,:),idx_static_1));
matches_static_13 = matches_13(:,ismember(matches_13(1,:),idx_static_1));
matches_static_14 = matches_14(:,ismember(matches_14(1,:),idx_static_1));

figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches_static_12, 'Stacking', 'v');
figure;
plotmatches(im1, im3, points_1(1:2,:), points_3(1:2,:), matches_static_13, 'Stacking', 'v');
figure;
plotmatches(im1, im4, points_1(1:2,:), points_4(1:2,:), matches_static_14, 'Stacking', 'v');

%% Compute fundamental matrices
p1 = [points_1(1:2, matches_static_12(1,:)); ones(1, length(matches_static_12))];
p2 = [points_2(1:2, matches_static_12(2,:)); ones(1, length(matches_static_12))];
[F2, ~] = ransac_fundamental_matrix(p1, p2, 2.0);

p1 = [points_1(1:2, matches_static_13(1,:)); ones(1, length(matches_static_13))];
p2 = [points_3(1:2, matches_static_13(2,:)); ones(1, length(matches_static_13))];
[F3, ~] = ransac_fundamental_matrix(p1, p2, 2.0);

p1 = [points_1(1:2, matches_static_14(1,:)); ones(1, length(matches_static_14))];
p2 = [points_4(1:2, matches_static_14(2,:)); ones(1, length(matches_static_14))];
[F4, ~] = ransac_fundamental_matrix(p1, p2, 2.0);

vgg_gui_F(im1rgb, im2rgb, F2');
vgg_gui_F(im1rgb, im3rgb, F3');
vgg_gui_F(im1rgb, im4rgb, F4');

%% Compute and show shared keypoints. Used to select the idx_I1
idx_shared = intersect(matches_12(1,:), matches_13(1,:));
idx_shared = intersect(idx_shared, matches_14(1,:));

figure;
imshow(im1rgb);
hold on;
plot(points_1(1,idx_shared), points_1(2,idx_shared), '+y');
text(points_1(1,idx_shared), points_1(2,idx_shared), compose('%d', idx_shared), 'Color', 'w');

%% Compute and plot trajectory

idx_I1 = 2022; % id obtained from idx_shared
idx_I2 = matches_12(2, (find(matches_12(1,:) == idx_I1)));
idx_I3 = matches_13(2, (find(matches_13(1,:) == idx_I1)));
idx_I4 = matches_14(2, (find(matches_14(1,:) == idx_I1)));

% Compute trajectory
point1_1 = [points_1(1:2, idx_I1)' 1]';
point1_2 = [points_2(1:2, idx_I2)' 1]';

l1 = cross(point1_1, point1_2); 
figure;imshow(im1rgb);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(points_1(1,idx_I1), points_1(2,idx_I1), 'y*');

% Compute the projection of point idx_I2 in the reference image
point2 = [points_2(1:2, idx_I2)' 1]';
l2 = F2'*point2;
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'c');
pi2 = cross(l1, l2);
plot(pi2(1)/pi2(3), pi2(2)/pi2(3), 'c*');

% Compute the projection of point idx_I3 in the reference image
point3 = [points_3(1:2, idx_I3)' 1]';
l3 = F3'*point3;
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
pi3 = cross(l1,l3);
plot(pi3(1)/pi3(3), pi3(2)/pi3(3), 'b*');

% Compute the projection of point idx_I4 in the reference image
point4 = [points_4(1:2,idx_I4)' 1]';
l4 = F4'*point4;
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
pi4 = cross(l1,l4);
plot(pi4(1)/pi4(3), pi4(2)/pi4(3), 'g*');

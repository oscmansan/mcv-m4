%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 2: Image mosaics

addpath('sift');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute image correspondences

%% Open images

images_llanes = {'Data/llanes/llanes_a.jpg',
                 'Data/llanes/llanes_b.jpg',
                 'Data/llanes/llanes_c.jpg'};

images_castle = {'Data/castle_int/0016_s.png',
                 'Data/castle_int/0015_s.png',
                 'Data/castle_int/0014_s.png'};

images_site13 = {'Data/aerial/site13/frame00000.png',
                 'Data/aerial/site13/frame00002.png',
                 'Data/aerial/site13/frame00003.png'};

images_site22 = {'Data/aerial/site22/frame_00001.tif',
                 'Data/aerial/site22/frame_00018.tif',
                 'Data/aerial/site22/frame_00030.tif'};

imargb = imread(images_llanes{1});
imbrgb = imread(images_llanes{2});
imcrgb = imread(images_llanes{3});

ima = double(rgb2gray(imargb)) / 255;
imb = double(rgb2gray(imbrgb)) / 255;
imc = double(rgb2gray(imcrgb)) / 255;

%% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01, 'Verbosity', 1);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01, 'Verbosity', 1);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01, 'Verbosity', 1);

fig1 = figure(1);
ax1 = axes('Parent', fig1);
imshow(imargb, [], 'Parent', ax1);
hold(ax1, 'on');
plot(ax1, points_a(1,:), points_a(2,:),'+y');
fig2 = figure(2);
ax2 = axes('Parent', fig2);
imshow(imbrgb, [], 'Parent', ax2);
hold(ax2, 'on');
plot(ax2, points_b(1,:), points_b(2,:),'+y');
fig3 = figure(3);
ax3 = axes('Parent', fig3);
imshow(imcrgb, [], 'Parent', ax3);
hold(ax3, 'on');
plot(ax3, points_c(1,:), points_c(2,:),'+y');

%% Match SIFT keypoints 

% between a and b
matches_ab = siftmatch(desc_a, desc_b);
figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab, 'Stacking', 'v');

% between b and c
matches_bc = siftmatch(desc_b, desc_c);
figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc, 'Stacking', 'v');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Compute the homography (DLT algorithm) between image pairs

%% Compute homography (normalized DLT) between a and b, play with the homography
th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function

figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab(:,inliers_ab), 'Stacking', 'v');

vgg_gui_H(imargb, imbrgb, Hab);

%% Compute homography (normalized DLT) between b and c, play with the homography
xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 

figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc(:,inliers_bc), 'Stacking', 'v');

vgg_gui_H(imbrgb, imcrgb, Hbc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Build the mosaic
corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, eye(3), corners);    % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab, corners);       % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);  % ToDo: complete the call to the function

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');

%% 3.1: ToDo: compute the mosaic with castle_int images
th = 3;
corners = [-600 1600 -250 800];
mosaic_castle = compute_mosaic_3(images_castle, th, corners);

figure;
imshow(mosaic_castle);
title('Castle mosaic A-B-C');

%% 3.2: ToDo: compute the mosaic with aerial images set 13
th = 3;
corners = [-300 1300 -100 1000];
mosaic_site13 = compute_mosaic_3(images_site13, th, corners);

figure;
imshow(mosaic_site13);
title('Site13 mosaic A-B-C');

%% 3.3: ToDo: compute the mosaic with aerial images set 22
th = 3;
corners = [-500 1500 -100 1100];
mosaic_site22 = compute_mosaic_3(images_site22, th, corners);

figure;
imshow(mosaic_site22);
title('Site22 mosaic A-B-C');

%% 3.4: ToDo: comment the results in every of the four cases: hypothetise why it works or
%       does not work

% COMMENTS PRESENTED IN THE LAB REPORT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Refine the homography with the Gold Standard algorithm

% Homography ab
x = points_a(1:2, matches_ab(1,inliers_ab));  % ToDo: set the non-homogeneous point coordinates of the 
xp = points_b(1:2, matches_ab(2,inliers_ab)); % point correspondences we will refine with the geometric method
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hab(:) ; x(:) ];      % The parameters or independent variables

Y_initial = gs_errfunction( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
err_initial = sum(sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'iter');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hab_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);

%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm
xhat = reshape(P(9+1:end), 2, []);
xhatp = euclid(Hab_r*[xhat; ones(1, length(xhat))]);

fig1 = figure(1);
ax1 = axes('Parent', fig1);
imshow(imargb, [], 'Parent', ax1);
hold(ax1, 'on');
plot(ax1, x(1,:), x(2,:),'+y');
plot(ax1, xhat(1,:), xhat(2,:),'+c');

fig2 = figure(2);
ax2 = axes('Parent', fig2);
imshow(imbrgb, [], 'Parent', ax2);
hold(ax2, 'on');
plot(ax2, xp(1,:), xp(2,:),'+y');
plot(ax2, xhatp(1,:), xhatp(2,:),'+c');

%%  Homography bc

% ToDo: refine the homography bc with the Gold Standard algorithm
x = points_b(1:2, matches_bc(1,inliers_bc));
xp = points_c(1:2, matches_bc(2,inliers_bc));
Xobs = [ x(:) ; xp(:) ];
P0 = [ Hbc(:) ; x(:) ];

options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'iter');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hbc_r = reshape( P(1:9), 3, 3 );

%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm
xhat = reshape(P(9+1:end), 2, []);
xhatp = euclid(Hbc_r*[xhat; ones(1, length(xhat))]);

fig1 = figure(1);
ax1 = axes('Parent', fig1);
imshow(imbrgb, [], 'Parent', ax1);
hold(ax1, 'on');
plot(ax1, x(1,:), x(2,:),'+y');
plot(ax1, xhat(1,:), xhat(2,:),'+c');

fig2 = figure(2);
ax2 = axes('Parent', fig2);
imshow(imcrgb, [], 'Parent', ax2);
hold(ax2, 'on');
plot(ax2, xp(1,:), xp(2,:),'+y');
plot(ax2, xhatp(1,:), xhatp(2,:),'+c');

%% Build mosaic
corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, eye(3), corners); % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab_r, corners); % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, inv(Hbc_r), corners); % ToDo: complete the call to the function

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Calibration with a planar pattern

clear all;

%% Read template and images.
T     = imread('Data/calib/template.jpg');
I{1}  = imread('Data/calib/graffiti1.tif');
I{2}  = imread('Data/calib/graffiti2.tif');
I{3}  = imread('Data/calib/graffiti3.tif');
%I{4}  = imread('Data/calib/graffiti4.tif');
%I{5}  = imread('Data/calib/graffiti5.tif');
Tg = sum(double(T), 3) / 3 / 255;
Ig{1} = sum(double(I{1}), 3) / 3 / 255;
Ig{2} = sum(double(I{2}), 3) / 3 / 255;
Ig{3} = sum(double(I{3}), 3) / 3 / 255;

N = length(I);

%% Compute keypoints.
fprintf('Computing sift points in template... ');
[pointsT, descrT] = sift(Tg, 'Threshold', 0.05, 'Verbosity', 1);
fprintf(' done\n');

points = cell(N,1);
descr = cell(N,1);
for i = 1:N
    fprintf('Computing sift points in image %d... ', i);
    [points{i}, descr{i}] = sift(Ig{i}, 'Threshold', 0.05, 'Verbosity', 1);
    fprintf(' done\n');
end

%% Match and compute homographies.
H = cell(N,1);
for i = 1:N
    % Match against template descriptors.
    fprintf('Matching image %d... ', i);
    matches = siftmatch(descrT, descr{i});
    fprintf('done\n');

    % Fit homography and remove outliers.
    x1 = [pointsT(1:2, matches(1,:)); ones(1, length(matches))];
    x2 = [points{i}(1:2, matches(2,:)); ones(1, length(matches))];
    H{i} = 0;
    [H{i}, inliers] =  ransac_homography_adaptive_loop(x1, x2, 3, 1000);

    % Plot inliers.
    figure;
    plotmatches(Tg, Ig{i}, pointsT(1:2,:), points{i}(1:2,:), matches(:, inliers));

    % Play with the homography
    %vgg_gui_H(T, I{i}, H{i});
end

%% Compute the Image of the Absolute Conic

Vij = @(H,i,j) [H(1,i)*H(1,j), H(1,i)*H(2,j)+H(2,i)*H(1,j), ...
                H(1,i)*H(3,j)+H(3,i)*H(1,j), H(2,i)*H(2,j), ...
                H(2,i)*H(3,j)+H(3,i)*H(2,j), H(3,i)*H(3,j)];
V = zeros(2*N, 6);
for i = 1:N
    V(2*i-1,:) = Vij(H{i},1,2);
    V(2*i,:) = Vij(H{i},1,1) - Vij(H{i},2,2);
end

[~,~,V] = svd(V);
o = V(:,end);
w = [o(1), o(2), o(3);
     o(2), o(4), o(5);
     o(3), o(5), o(6)];
 
%% Recover the camera calibration.

KKt = inv(w);
K = chol(KKt);  % shouldn't we use the lower triangular?
    
% ToDo: in the report make some comments related to the obtained internal
%       camera parameters and also comment their relation to the image size

% COMMENTS PRESENTED IN THE LAB REPORT

%% Compute camera position and orientation.
R = cell(N,1);
t = cell(N,1);
P = cell(N,1);
figure;hold;
for i = 1:N
    % ToDo: compute r1, r2, and t{i}
    r1 = K \ H{i}(:,1);
    r2 = K \ H{i}(:,2);
    t{i} = K \ H{i}(:,3);
    
    % Solve the scale ambiguity by forcing r1 and r2 to be unit vectors.
    s = sqrt(norm(r1) * norm(r2)) * sign(t{i}(3));
    r1 = r1 / s;
    r2 = r2 / s;
    t{i} = t{i} / s;
    R{i} = [r1, r2, cross(r1,r2)];
    
    % Ensure R is a rotation matrix
    [U, S, V] = svd(R{i});
    R{i} = U * eye(3) * V';
   
    P{i} = K * [R{i}, t{i}];
    plot_camera(P{i}, 800, 600, 200);
end

% ToDo: in the report explain how the optical center is computed in the
%       provided code

% COMMENTS PRESENTED IN THE LAB REPORT

[ny,nx] = size(T);
p1 = [-nx/2 -ny/2 0]';
p2 = [nx/2 -ny/2 0]';
p3 = [nx/2 ny/2 0]';
p4 = [-nx/2 ny/2 0]';
% Draw planar pattern
vgg_scatter_plot([p1 p2 p3 p4 p1], 'g');
% Paint image texture
surface('XData',[-nx/2, nx/2; -nx/2, nx/2],'YData',[0 0; 0 0],'ZData',[ny/2, ny/2; -ny/2, -ny/2],'CData',T,'FaceColor','texturemap');
colormap(gray);
axis equal;

%% Plot a static camera with moving calibration pattern.
figure; hold;
plot_camera(K * eye(3,4), 800, 600, 200);
% ToDo: complete the call to the following function with the proper
%       coordinates of the image corners in the new reference system
corners = [p1, p2, p3, p4, p1; ones(1,5)];
for i = 1:N
    vgg_scatter_plot([R{i},t{i}] * corners, 'r');
end

%% Augmented reality: Plot some 3D points on every camera.
[Th, Tw] = size(Tg);
cube = [0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 0 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 0 0; 1 0 1; 1 0 1; 0 0 1; 0 0 1; 0 0 0; 0 1 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 0; 0 0 0; 0 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 1 0 1; 1 0 1; 1 0 0 ]';

X = (cube - .5) * Tw / 4 + repmat([Tw / 2; Th / 2; -Tw / 8], 1, length(cube));
X = [X; ones(1, length(X))];

for i = 1:N
    figure; colormap(gray);
    imagesc(Ig{i});
    hold on;
    x = euclid(P{i} * X);
    vgg_scatter_plot(x, 'g');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. OPTIONAL: Detect the UPF logo in the two UPF images using the 
%%              DLT algorithm (folder "logos").
%%              Interpret and comment the results.

% Manual Keypoints Detection
img_dst = imread('Data/logos/UPFbuilding.jpg');
img_src = imread('Data/logos/logo_master.png');

% corners coordinates of logo in the building image
pts_dst = [321, 423, 310, 416;
           47,   65, 123, 134]; 

[h,w,c] = size(img_src);
pts_src = [0, w, w,  0;
           0,  0,   h, h];

figure;
imshow(img_dst);
hold on;
plot(pts_dst(1,:), pts_dst(2,:),'+y');
hold off;

%TODO: complete with match from manual keypoints detection

% Auto Keypoints Detection
img_dst_auto = imread('Data/logos/UPFstand.jpg');
img_match_auto = imread('Data/logos/logoUPF.png');
img_src_auto = imread('Data/logos/logo_master.png');

[rows, cols, a] = size(img_match_auto);
img_src_auto =  imresize(img_src_auto, [rows cols]);

img_dest_gray = double(rgb2gray(img_dst_auto)) / 255;
img_match_gray = double(rgb2gray(img_match_auto)) / 255;

% Compute and match SIFT keypoints
[points_a, desc_a] = sift(img_match_gray, 'Threshold', 0.01);
[points_b, desc_b] = sift(img_dest_gray, 'Threshold', 0.01);
matches_ab = siftmatch(desc_a, desc_b);

%Fit homography and remove outliers
pts_src_auto = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
pts_dst_auto = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(pts_src_auto, pts_dst_auto, 3, 1000);

figure;
plotmatches(img_match_gray, img_dest_gray, points_a(1:2,:), points_b(1:2,:), matches_ab(:,inliers_ab), 'Stacking', 'v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. OPTIONAL: Replace the logo of the UPF by the master logo
%%              in one of the previous images using the DLT algorithm.

% Automatic Replacement
[h,w,c] = size(img_dst_auto);
corners = [0 w-1 0 h-1];

warp_src = apply_H_v2(img_src_auto, Hab, corners);
merge = max(img_dst_auto, warp_src);
figure;
imshow(merge)


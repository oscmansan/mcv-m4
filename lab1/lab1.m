%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
I=imread('Data/0005_s.png'); % we have to be in the proper folder
I = imresize(I,0.5);

% ToDo: generate a matrix H which produces a similarity transformation

theta = 5/180*pi;
scale = 0.5;
trans = [10,20];
H = [scale*cos(theta),  scale*-sin(theta),  trans(1);
     scale*sin(theta),  scale*cos(theta),  trans(2);
     0,  0,  1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation

% PROVA INICIAL
% vert_angle_rotation = 50/180*pi;
% hor_angle_rotation = 5/180*pi;
% vertical_scaling = 1;
% horizontal_scaling = 0.5;
% translation_vector = [10,20];
% H = [vertical_scaling*cos(vert_angle_rotation),  horizontal_scaling*-sin(hor_angle_rotation),  translation_vector(1);
%      vertical_scaling*sin(vert_angle_rotation),  horizontal_scaling*cos(hor_angle_rotation),  translation_vector(2);
%      0,  0,  1];

theta = 30/180*pi;
phi = 15/180*pi;
scale = 0.5;
trans = [10,20];
H = [scale*(cos(theta)*cos(phi)-sin(theta)*sin(phi)),  scale*(-cos(theta)*sin(phi)-sin(theta)*cos(phi)),  scale*trans(1)*(cos(theta)*cos(phi)-sin(theta)*sin(phi))+scale*trans(2)*(-cos(theta)*sin(phi)-sin(theta)*cos(phi));
     scale*(sin(theta)*cos(phi)+cos(theta)*sin(phi)),  scale*(-sin(theta)*sin(phi)+cos(theta)*cos(phi)),  scale*trans(1)*(sin(theta)*cos(phi)+cos(theta)*sin(phi))+scale*trans(2)*(-sin(theta)*sin(phi)+cos(theta)*cos(phi));
     0,  0,  1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation

Hrot1 = [cos(theta),  -sin(theta),  0;
         sin(theta),  cos(theta),  0;
         0,  0,  1];

Hrot2 = [cos(phi),  -sin(phi),  0;
         sin(phi),  cos(phi),  0;
         0,  0,  1];

Hscale = [scale,  0,  0;
          0,  scale,  0;
          0,  0,  1];

Htrans = [1,  0,  trans(1);
          0,  1,  trans(2);
          0,  0,  1];

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above

if sum(sum(Hrot1*Hrot2*Hscale*Htrans-H)) < 1e-10 % Can't correctly compare floating points directly
    display('CORRECT MATRIX DECOMPOSITION')
else
    display('WRONG MATRIX DECOMPOSITION')
end

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before

I3 = apply_H(I, Hrot1);
I3 = apply_H(I3, Hrot2);
I3 = apply_H(I3, Hscale);
I3 = apply_H(I3, Htrans);
figure; imshow(uint8(I3));

% TODO: Aqu� haur�em de comparar num�ricament les dues imatges

%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation
Hproj = H;
Hproj(3,:)=[0, 0.003, 1];

I2 = apply_H(I, Hproj);
figure; imshow(I); figure; imshow(uint8(I2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification
% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = get_line_from_points(p1,p2);
l2 = get_line_from_points(p3,p4);
l3 = get_line_from_points(p5,p6);
l4 = get_line_from_points(p7,p8);


% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image

I2 = apply_H(I, H);
figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that are used in every step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)




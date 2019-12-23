close all;
clear;
clc;
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

% ToDo: generate a matrix H which produces a similarity transformation
theta = 30/180*pi;
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
s = 0.5;
t = [10,20];
H = [s*R(1,1),  s*R(1,2),  t(1);
     s*R(2,1),  s*R(2,2),  t(2);
     0,  0,  1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation
A = [0.5, -0.25; 0.25, 0.5];
t = [10,20];
H = [A(1,1),  A(1,2),  t(1);
     A(2,1),  A(2,2),  t(2);
     0,  0,  1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
[U,D,V] = svd(A);  % A = (U*V')*(V*D*V') = R1*R2'*D*R2
R1 = U*V';
R2 = V';

Hrot1 = [R1(1,1), R1(1,2), 0;
         R1(2,1), R1(2,2), 0;
         0, 0, 1];
     
Hrot2 = [R2(1,1), R2(1,2), 0;
         R2(2,1), R2(2,2), 0;
         0, 0, 1];
     
Hscale = [D(1,1), D(1,2), 0;
          D(2,1), D(2,2), 0;
          0, 0, 1];
      
Htrans = [1, 0, t(1);
          0, 1, t(2);
          0, 0, 1];

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above
H2 = Htrans*Hrot1*Hrot2'*Hscale*Hrot2;
assert(max(abs(H2-H),[],'all') < 1e-10);

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
I3 = apply_H(I, Hrot2);
I3 = apply_H(I3, Hscale);
I3 = apply_H(I3, Hrot2');
I3 = apply_H(I3, Hrot1);
I3 = apply_H(I3, Htrans);
figure; imshow(uint8(I3));  % the result is blurrier due to successive interpolations


%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation
A = [0.5, -0.25; 0.25, 0.5];
t = [10,20];
v = [0,0.0005];
H = [A(1,1),  A(1,2),  t(1);
     A(2,1),  A(2,2),  t(2);
     v(1),  v(2),  1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification
% choose the image points
I = imread('Data/0000_s.png');
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
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5,p6);
l4 = cross(p7,p8);
l5 = cross(p1,p4);
l6 = cross(p2,p3);

% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');
hold off;

% ToDo: compute the homography that affinely rectifies the image
v1 = cross(l1,l2);
v2 = cross(l3,l4);
linf = cross(v1,v2);
linf = linf / linf(3);  % why does this need to be normalized?
 
Hpa = [1, 0, 0; 0, 1, 0; linf'];
 
I2 = apply_H(I, Hpa);
figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
lr1 = Hpa'\l1;
lr2 = Hpa'\l2;
lr3 = Hpa'\l3;
lr4 = Hpa'\l4;
lr5 = Hpa'\l5;
lr6 = Hpa'\l6;

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'g');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'g');
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'c');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'c');
hold off;

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation
angle_l1_l2 = angle_between_lines(l1,l2)/pi*180;
angle_lr1_lr2 = angle_between_lines(lr1,lr2)/pi*180;

angle_l3_l4 = angle_between_lines(l3,l4)/pi*180;
angle_lr3_lr4 = angle_between_lines(lr3,lr4)/pi*180;

angle_l1_l3 = angle_between_lines(l1,l3)/pi*180;
angle_l2_l4 = angle_between_lines(l2,l4)/pi*180;

fprintf('Angle between (l1, l2): %.2f -> %.2f\n', angle_l1_l2, angle_lr1_lr2);
fprintf('Angle between (l3, l4): %.2f -> %.2f\n', angle_l3_l4, angle_lr3_lr4);

assert(angle_lr1_lr2 == 0);
assert(angle_lr3_lr4 == 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

lm1 = [lr1(1)*lr3(1), lr1(1)*lr3(2) + lr1(2)*lr3(1), lr1(2)*lr3(2)];
%lm2 = [lr2(1)*lr4(1), lr2(1)*lr4(2) + lr2(2)*lr4(1), lr2(2)*lr4(2)];
lm2 = [lr5(1)*lr6(1), lr5(1)*lr6(2) + lr5(2)*lr6(1), lr5(2)*lr6(2)];

A = [lm1; lm2];  % need to have different directions
s = null(A);

S = [s(1), s(2); s(2), s(3)];
K = chol(S, 'lower');

Has = [inv(K) [0;0]; 0, 0, 1];
H3 = Has*Hpa;
I3 = apply_H(I, H3);

% compute the lines from the transformed points
lrr1 = H3'\l1;
lrr2 = H3'\l2;
lrr3 = H3'\l3;
lrr4 = H3'\l4;
lrr5 = H3'\l5;
lrr6 = H3'\l6;

figure; imshow(uint8(I3));
hold on;
t=1:0.1:1000;
plot(t, -(lrr1(1)*t + lrr1(3)) / lrr1(2), 'y');
plot(t, -(lrr2(1)*t + lrr2(3)) / lrr2(2), 'y');
plot(t, -(lrr3(1)*t + lrr3(3)) / lrr3(2), 'g');
plot(t, -(lrr4(1)*t + lrr4(3)) / lrr4(2), 'g');
plot(t, -(lrr5(1)*t + lrr5(3)) / lrr5(2), 'c');
plot(t, -(lrr6(1)*t + lrr6(3)) / lrr6(2), 'c');
hold off;

% compute the angle between pairs of lines after rectification
angle_lr1_lr3 = angle_between_lines(lrr1,lrr3)/pi*180;
angle_lr2_lr4 = angle_between_lines(lrr2,lrr4)/pi*180;
angle_lr5_lr6 = angle_between_lines(lrr5,lrr6)/pi*180;

fprintf('Angle between (lr1, lr3): %.2f\n', angle_lr1_lr3);
fprintf('Angle between (lr2, lr4): %.2f\n', angle_lr2_lr4);
fprintf('Angle between (lr5, lr6): %.2f\n', angle_lr5_lr6);

assert(angle_lr1_lr3 == 90);
assert(angle_lr2_lr4 == 90);
assert(angle_lr5_lr6 == 90);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that are used in every step.

% choose the image points
I = imread('Data/0001_s.png');
A = load('Data/0001_s_info_lines.txt');

% crop the image
I = I(:, 1:450, :);

% indices of lines
i = 614;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 159;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 645;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 541;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1,p2);
l2 = cross(p3,p4);
l3 = cross(p5,p6);
l4 = cross(p7,p8);

% crossed lines
p9 = cross(l1, l3);
p10 = cross(l1, l4);
p11 = cross(l2, l3);
p12 = cross(l2, l4);
l5 = cross(p9, p12);
l6 = cross(p10, p11);

% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'y');
plot(t, -(l6(1)*t + l6(3)) / l6(2), 'y');
hold off;

% compute the homography that affinely rectifies the image
v1 = cross(l1,l2);
v2 = cross(l3,l4);
 
linf = cross(v1, v2);
linf = linf/linf(3);

Hpa = [1 0 0; 0 1 0; linf'];
 
I2 = apply_H(I, Hpa);
figure; imshow(uint8(I2));

% compute the transformed lines lr1, lr2, lr3, lr4
lr1 = Hpa'\l1;
lr2 = Hpa'\l2;
lr3 = Hpa'\l3;
lr4 = Hpa'\l4;
lr5 = Hpa'\l5;
lr6 = Hpa'\l6;

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'g');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'g');
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'c');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'c');
hold off;

% compute the homography that does metric rectification
lm1 = [lr1(1)*lr3(1), lr1(1)*lr3(2) + lr1(2)*lr3(1), lr1(2)*lr3(2)];
%lm2 = [lr2(1)*lr4(1), lr2(1)*lr4(2) + lr2(2)*lr4(1), lr2(2)*lr4(2)];
lm2 = [lr5(1)*lr6(1), lr5(1)*lr6(2) + lr5(2)*lr6(1), lr5(2)*lr6(2)];

A = [lm1; lm2];
s = null(A);

S = [s(1), s(2); s(2), s(3)];
K = chol(S, 'lower');

Has = [inv(K) [0;0]; 0 0 1];
[I3, minX, minY] = apply_H(I2, Has);

% compute the lines from the transformed points
lr1 = Has'\lr1;
lr2 = Has'\lr2;
lr3 = Has'\lr3;
lr4 = Has'\lr4;
lr5 = Has'\lr5;
lr6 = Has'\lr6;

figure;imshow(uint8(I3));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)+(minY-1)*lr1(2)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)+(minY-1)*lr2(2)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)+(minY-1)*lr3(2)) / lr3(2), 'g');
plot(t, -(lr4(1)*t + lr4(3)+(minY-1)*lr4(2)) / lr4(2), 'g');
plot(t, -(lr5(1)*t + lr5(3)+(minY-1)*lr5(2)) / lr5(2), 'c');
plot(t, -(lr6(1)*t + lr6(3)+(minY-1)*lr6(2)) / lr6(2), 'c');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)

I = imread('Data/0001_s.png');
A = load('Data/0001_s_info_lines.txt');

I = I(:, 1:450, :);

i = 315;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
l1 = cross(p1,p2);
i = 301;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
m1 = cross(p1,p2);

i = 614;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
l2 = cross(p1,p2);
i = 541;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
m2 = cross(p1,p2);

i = 188;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
l3 = cross(p1,p2);
i = 337;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
m3 = cross(p1,p2);

i = 284;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
l4 = cross(p1,p2);
i = 298;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
m4 = cross(p1,p2);

i = 645;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
l5 = cross(p1,p2);
i = 170;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
m5 = cross(p1,p2);

figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(m1(1)*t + m1(3)) / m1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'g');
plot(t, -(m2(1)*t + m2(3)) / m2(2), 'g');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'c');
plot(t, -(m3(1)*t + m3(3)) / m3(2), 'c');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'b');
plot(t, -(m4(1)*t + m4(3)) / m4(2), 'b');
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'r');
plot(t, -(m5(1)*t + m5(3)) / m5(2), 'r');
hold off;

%waitforbuttonpress;

lm1 = [l1(1)*m1(1), (l1(1)*m1(2)+l1(2)*m1(1))/2, l1(2)*m1(2), (l1(1)*m1(3)+l1(3)*m1(1))/2, (l1(2)*m1(3)+l1(3)*m1(2))/2, l1(3)*m1(3)];
lm2 = [l2(1)*m2(1), (l2(1)*m2(2)+l2(2)*m2(1))/2, l2(2)*m2(2), (l2(1)*m2(3)+l2(3)*m2(1))/2, (l2(2)*m2(3)+l2(3)*m2(2))/2, l2(3)*m2(3)];
lm3 = [l3(1)*m3(1), (l3(1)*m3(2)+l3(2)*m3(1))/2, l3(2)*m3(2), (l3(1)*m3(3)+l3(3)*m3(1))/2, (l3(2)*m3(3)+l3(3)*m3(2))/2, l3(3)*m3(3)];
lm4 = [l4(1)*m4(1), (l4(1)*m4(2)+l4(2)*m4(1))/2, l4(2)*m4(2), (l4(1)*m4(3)+l4(3)*m4(1))/2, (l4(2)*m4(3)+l4(3)*m4(2))/2, l4(3)*m4(3)];
lm5 = [l5(1)*m5(1), (l5(1)*m5(2)+l5(2)*m5(1))/2, l5(2)*m5(2), (l5(1)*m5(3)+l5(3)*m5(1))/2, (l5(2)*m5(3)+l5(3)*m5(2))/2, l5(3)*m5(3)];

A = [lm1; lm2; lm3; lm4; lm5];
s = null(A);  % should return a single solution
 
a = s(1);
b = s(2);
c = s(3);
d = s(4);
e = s(5);
f = s(6);
C = [a, b/2, d/2; b/2, c, e/2; d/2, e/2, f];

assert(issymmetric(C));

[U,D,V] = svd(C);
H = U;

lr1 = H'\l1;
mr1 = H'\m1;
lr2 = H'\l2;
mr2 = H'\m2;
lr3 = H'\l3;
mr3 = H'\m3;
lr4 = H'\l4;
mr4 = H'\m4;
lr5 = H'\l5;
mr5 = H'\m5;

[I2, minX, minY] = apply_H(I, H);

figure;imshow(uint8(I2));
% hold on;
% t=1:0.1:1000;
% plot(t, -(lr1(1)*t + lr1(3) +(minY-1)*lr1(2)) / lr1(2), 'y');
% plot(t, -(mr1(1)*t + mr1(3) +(minY-1)*mr1(2)) / mr1(2), 'y');
% plot(t, -(lr2(1)*t + lr2(3) +(minY-1)*lr2(2)) / lr2(2), 'g');
% plot(t, -(mr2(1)*t + mr2(3) +(minY-1)*mr2(2)) / mr2(2), 'g');
% plot(t, -(lr3(1)*t + lr3(3) +(minY-1)*lr3(2)) / lr3(2), 'c');
% plot(t, -(mr3(1)*t + mr3(3) +(minY-1)*mr3(2)) / mr3(2), 'c');
% plot(t, -(lr4(1)*t + lr4(3) +(minY-1)*lr4(2)) / lr4(2), 'b');
% plot(t, -(mr4(1)*t + mr4(3) +(minY-1)*mr4(2)) / mr4(2), 'b');
% plot(t, -(lr5(1)*t + lr5(3) +(minY-1)*lr5(2)) / lr5(2), 'r');
% plot(t, -(mr5(1)*t + mr5(3) +(minY-1)*mr5(2)) / mr5(2), 'r');
% hold off;

fprintf('angle between (l1,m1): %.3f\n', angle_between_lines(lr1,mr1)/pi*180);
fprintf('angle between (l2,m2): %.3f\n', angle_between_lines(lr2,mr2)/pi*180);
fprintf('angle between (l3,m3): %.3f\n', angle_between_lines(lr3,mr3)/pi*180);
fprintf('angle between (l4,m4): %.3f\n', angle_between_lines(lr4,mr4)/pi*180);
fprintf('angle between (l5,m5): %.3f\n', angle_between_lines(lr5,mr5)/pi*180);

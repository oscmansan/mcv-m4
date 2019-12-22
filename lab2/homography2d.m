function [H] = homography2d(x1_points, x2_points)

% Normalise each set of points so that the origin 
% is at centroid and mean distance from origin is sqrt(2).
% https://es.mathworks.com/matlabcentral/fileexchange/54544-normalise2dpts-pts
[x1_norm, T1] = normalise2dpts(x1_points);
[x2_norm, T2] = normalise2dpts(x2_points);

% Assemble the 2nx9 matrices A_i
n_points = length(x1_norm);
A = zeros(2*n_points,9);

for n = 1:n_points
    X = x1_norm(:,n)';
    x = x2_norm(1,n); 
    y = x2_norm(2,n); 
    w = x2_norm(3,n);
    A(2*n-1,:) = [[0, 0, 0], -w*X, y*X];
    A(2*n,:) = [w*X, [0, 0, 0], -x*X];
end

% Singular Value Decomposition
% http://es.mathworks.com/help/matlab/ref/svd.html?s_tid=srchtitle
[~,~,V] =svd(A);

% Take the last column of the transposed of V, that's the singular
% vector with the lowest singular value.
h = V(:,9);

% Reshape h to be a 3x3 matrix.
H = reshape(h,3,3)';

% Denormalization
H = T2\H*T1;

return
function F = fundamental_matrix(x1, x2)

assert(length(x1) == length(x2));
N = length(x1);

% Normalize: compute H and H'
[x1, H1] = normalise2dpts(x1);
[x2, H2] = normalise2dpts(x2);

% Apply 8-point algorithm to estimate F

% Create matrix W from points correspondences
W = zeros(N,9);
for i=1:N
    W(i,:) = [x1(1,i)*x2(1,i), x1(2,i)*x2(1,i), x2(1,i),...
              x1(1,i)*x2(2,i), x1(2,i)*x2(2,i), x2(2,i),...
              x1(1,i), x1(2,i), x1(3,i)];
end

% Compute SVD of W
[~, ~, V] = svd(W);
% Create vector f from last column of V
f = V(:,end);
% Compose fundamental matrix F_rank3
F_rank3 = reshape(f,[3,3]).';
% Remove last singular value of D to create ^D
[U,D,V] = svd(F_rank3);
D(end,end) = 0;
% Recompute matrix F_rank2
F_rank2 = U * D * V';

% Denormalize F
F = H2.' * F_rank2 * H1;

end
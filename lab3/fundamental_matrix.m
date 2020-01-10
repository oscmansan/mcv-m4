function F = fundamental_matrix(x1, x2)

% Normalize: compute H and H'
[x1_norm, H1] = normalise2dpts(x1);
[x2_norm, H2] = normalise2dpts(x2);

% Apply 8-point algorithm to estimate F

    % Create matrix W from points correspondences
    for i=1:size(x1_norm,2)
        W(i,:) = [x1_norm(1,i)*x2_norm(1,i), x1_norm(2,i)*x2_norm(1,i), x2_norm(1,i),...
                  x1_norm(1,i)*x2_norm(2,i), x1_norm(2,i)*x2_norm(2,i), x2_norm(2,i),...
                  x1_norm(1,i), x1_norm(2,i), x1_norm(3,i)];
    end
    % Compute SVD of W
    [~, ~, V] = svd(W);
    % Create vector f from last column of V
    f = V(:,end);
    % Compose fundamental matrix F_rank3
    F_rank3 = [f(1:3),f(4:6),f(7:9)]';
    % Remove last singular value of D to create ^D
    [U,D,V] = svd(F_rank3);
    D(end,end) = 0;
    % Recompute matrix F_rank2
    F_rank2 = U * D * V';

% Denormalize F
F = H2' * F_rank2 * H1;

end
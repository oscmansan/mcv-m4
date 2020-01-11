function [F, idx_inliers] = ransac_fundamental_matrix(x1, x2, th, max_it)

if nargin < 4
    max_it = 1000;
end

[Ncoords, Npoints] = size(x1);

% ransac
it = 0;
best_inliers = [];
% probability that at least one random sample set is free of outliers
p = 0.999; 
while it < max_it    
    points = randsample(Npoints, 8);
    F = fundamental_matrix(x1(:,points), x2(:,points));
    inliers = compute_inliers(F, x1, x2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(inliers)/Npoints;
    pNoOutliers = 1 -  fracinliers^4;
    pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
    max_it = log(1-p)/log(pNoOutliers);
    
    it = it + 1;
end

% compute F from all the inliers
F = fundamental_matrix(x1(:,best_inliers), x2(:,best_inliers));
idx_inliers = best_inliers;

end

function idx_inliers = compute_inliers(F, x1, x2, th)

N = length(x1);

x2tFx1 = zeros(1,N);
for i = 1:N
    x2tFx1(i) = x2(:,i)'*F*x1(:,i);
end

Fx1 = F*x1;
Ftx2 = F.'*x2;

Fx1 = normalise(Fx1);
Ftx2 = normalise(Ftx2);

% Sampson distance
d2 = x2tFx1.^2 ./ (Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2);

idx_inliers = find(d2 < th^2);

end

function xn = normalise(x)    
    xn = x ./ repmat(x(end,:), size(x,1), 1);
end

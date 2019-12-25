function E = gs_errfunction(P, Xobs)
% GS_ERRFUNCTION Compute the Gold Standard error.
%   Arguments:
%       - P: [ H(:) ; x(:) ]
%       - Xobs: [ x(:) ; xp(:) ]
%   All points are received in inhomogeneous coordinates (2 dimensions)

% get estimated matrix H contained in P
H = reshape(P(1:9), 3, 3);

% get observed points x and xp contained in Xobs
nPoints = length(Xobs)/2;
x = reshape(Xobs(1:nPoints), 2, []);
xp = reshape(Xobs(nPoints+1:end), 2, []);

% get estimated points xhat contained in P and compute xhatp
xhat = reshape(P(9+1:end), 2, []);
xhatp = euclid(H*[xhat; ones(1, length(xhat))]);

% compute reprojection error
e1 = l2_dist(x, xhat);
e2 = l2_dist(xp, xhatp);
E = sqrt(e1+e2);  % because lsqnonlin implicitly squares each element

end

function d = l2_dist(a, b)
    d = sum((a-b).^2);
end
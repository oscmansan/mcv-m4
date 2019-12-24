function E = gs_errfunction(P, Xobs)

H = reshape(P(1:9), 3, 3);

nPoints = length(Xobs)/2;
x = reshape(Xobs(1:nPoints), 2, []);
xp = reshape(Xobs(nPoints+1:end), 2, []);

xhat = reshape(P(9+1:end), 2, []);
xhatp = H*[xhat; ones(1, length(xhat))];

e1 = l2_dist(x, euclid(xhat));
e2 = l2_dist(xp, euclid(xhatp));

E = e1+e2;

end

function d = l2_dist(x, y)
    d = (sum((x-y).^2));
end
function [error, mean_err] = reprojection_error(P1, P2, X, x1, x2)

    x1_reprojected = euclid(P1*X);
    x2_reprojected = euclid(P2*X);
    
    error = distance(x1, x1_reprojected) + distance(x2, x2_reprojected);

    mean_err = mean(error);
end

function [dist] = distance(x,y)
    dist = sqrt((x(1,:) - y(1,:)).^2 + (x(2,:) - y(2,:)).^2);
end
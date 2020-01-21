function [error, mean_err] = reprojection_error(P1, P2, X, x1, x2)

    x1_reprojected = euclid(P1*X);
    x2_reprojected = euclid(P2*X);
    
    error1 = distance(x1, x1_reprojected);
    error2 = distance(x2, x2_reprojected);
    
    error = error1 + error2;

    mean_err = (mean(error1)+mean(error2))/2;
    
    figure;
    histogram(error1)
    hold on
    histogram(error2)
    line([mean_err mean_err], ylim, 'Color','g');
    legend('hist error camera 1', 'hist error camera 2', 'mean error value');
    hold off        
    
end

function [dist] = distance(x,y)
    dist = sqrt(sum(x-y).^2);
end
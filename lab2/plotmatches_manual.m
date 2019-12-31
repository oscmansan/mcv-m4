function plotmatches_manual(im1, im2, points1, points2, inliers)
    [rows1,cols1] = size(im1);
    [rows2,cols2] = size(im2);

    % Create blank image
    stackedImage = zeros(max([rows1,rows2]), cols1+cols2);
    stackedImage = cast(stackedImage, class(im1)); %// Make sure we cast output
    % Place two images side by side
    stackedImage(1:rows1,1:cols1) = im1;
    stackedImage(1:rows2,cols1+1:cols1+cols2) = im2;

    % Code from before
    imshow(stackedImage);
    width = size(im1, 2);
    hold on;
    
    numPoints = length(inliers);
    for n = 1 : numPoints
        i = inliers(n);
        points1_x = points1(1, i);
        points1_y = points1(2, i);

        points2_x = points2(1,i);
        points2_y = points2(2, i);
        plot(points1_x,points1_y, 'g+',points2_x + width, ...
            points2_y, 'g+');
        line([points1_x points2_x + width], [points1_y points2_y], ...
             'Color', 'green');
    end
    
end

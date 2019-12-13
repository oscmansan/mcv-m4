function [I2] = apply_H(I,H)

I = double(I);
[height, width, nChannels] = size(I);

% TODO: Calcular la mida del resultat (Aquí està posat dummy)
height2 = height;
width2 = width;

% Initialize final image
I2 = zeros(height2,width2,nChannels);

% Define initial set of points
[X,Y] = meshgrid(1:width,1:height); % MATLAB gira la X i la Y sembla ser

% Initialize final set of points
[X2,Y2] = meshgrid(1:width2,1:height2);

% Compute final set of points
for i = 1:height2
    for j = 1:width2
        p2_point = [i;j;1];
        new_p2_point = H\p2_point; % Inverse matrix multiplication
        new_r2_point = new_p2_point(1:2)/new_p2_point(3);
        if (new_r2_point(1) < 0 || new_r2_point(1) > height || new_r2_point(2) < 0 || new_r2_point(2) > width)
            X2(i,j) = 0;
            Y2(i,j) = 0;
        else
            X2(i,j) = new_r2_point(2);
            Y2(i,j) = new_r2_point(1);
        end
    end
end

% Assign values to final pixels
for nChannel=1:nChannels
    I2(:,:,nChannel)=interp2(X,Y,I(:,:,nChannel),X2,Y2);
end

end
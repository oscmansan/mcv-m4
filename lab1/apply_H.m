function I2 = apply_H(I,H)

I = double(I);
[numRows, numCols, numChannels] = size(I);

% compute transformed image corners in world coordinates
topLeftH = H * [1; 1; 1];
topRightH = H * [numCols; 1; 1];
botRightH = H * [numCols; numRows; 1];
botLeftH = H * [1; numRows; 1];
cornersH = [topLeftH, topRightH, botRightH, botLeftH];
corners = [cornersH(1,:)./cornersH(3,:); cornersH(2,:)./cornersH(3,:)];

% compute transformed image limits in world coordinates
minX = ceil(min(corners(1,:)));
maxX = floor(max(corners(1,:)));
minY = ceil(min(cornersH(2,:)));
maxY = floor(max(corners(2,:)));

dstWidth = maxX-minX+1;
dstHeight = maxY-minY+1;

% compute mapping from dst pixels in world coordinates to src pixels
[dstX, dstY] = meshgrid(minX:maxX, minY:maxY);
dstPointsH = [dstX(:).'; dstY(:).'; ones(1, dstWidth*dstHeight)];
%srcPoints = inv(H) * dstPointsH;
srcPointsH = H\dstPointsH;  % faster
srcX = reshape(srcPointsH(1,:)./srcPointsH(3,:), [dstHeight,dstWidth]);
srcY = reshape(srcPointsH(2,:)./srcPointsH(3,:), [dstHeight,dstWidth]);

% copy values from src pixels to dst pixels in intrinsic coordinates
I2 = zeros(dstHeight, dstWidth, numChannels);
for k = 1:numChannels
    I2(:,:,k) = interp2(I(:,:,k), srcX, srcY);
end

end
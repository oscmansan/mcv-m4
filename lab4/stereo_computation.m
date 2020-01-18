function disparity = stereo_computation(Ileft,Iright,minDisp,maxDisp,winSize,cost)
%STEREO_COMPUTATION Summary of this function goes here
%   Detailed explanation goes here

[numRows,numCols] = size(Ileft);
p = floor(winSize/2);

switch cost
    case 'SSD'
        cost = @(x,y) norm(x(:)-y(:),2).^2;
    otherwise
        warning('Unexpected cost type.');
end

disparity = zeros(size(Ileft));
tic
for i = 1+p:numRows-p
    for j = 1+p:numCols-p
        cropLeft = Ileft(i-p:i+p,j-p:j+p);
        minCost = Inf;
        for k = [max(1+p,j-maxDisp):max(1+p,j-minDisp) min(numCols-p,j+minDisp):min(numCols-p,j+maxDisp)]
            cropRight = Iright(i-p:i+p,k-p:k+p);
            c = cost(cropLeft,cropRight);
            if c < minCost
                minCost = c;
                idx = k;
            end
        end
        disparity(i,j) = abs(j-idx);
    end
    fprintf('%.2f %%\n',(i-p)/(numRows-2*p)*100);
end
toc

end


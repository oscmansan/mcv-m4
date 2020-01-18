function disparity = stereo_computation(Ileft,Iright,minDisp,maxDisp,winSize,cost)
%STEREO_COMPUTATION Summary of this function goes here
%   Detailed explanation goes here

[numRows,numCols] = size(Ileft);
p = floor(winSize/2);

switch cost
    case 'SAD'
        cost = @sad;
    case 'SSD'
        cost = @ssd;
    case 'NCC'
        cost = @(x,y,w) 1-ncc(x,y,w);
    otherwise
        warning('Unexpected cost type.');
end

disparity = zeros(size(Ileft));
tic
for i = 1+p:numRows-p
    for j = 1+p:numCols-p
        cropLeft = Ileft(i-p:i+p,j-p:j+p);
        w = ones(winSize)/winSize^2;
        minCost = Inf;
        for k = [max(1+p,j-maxDisp):max(1+p,j-minDisp) min(numCols-p,j+minDisp):min(numCols-p,j+maxDisp)]
            cropRight = Iright(i-p:i+p,k-p:k+p);
            c = cost(cropLeft,cropRight,w);
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

function c = sad(x,y,w)
    c = sum(w(:).*abs(x(:)-y(:)));
end

function c = ssd(x,y,w)
    c = sum(w(:).*(x(:)-y(:)).^2);
end

function c = ncc(x,y,w)
    sumx = sum(w(:).*x(:));
    sumy = sum(w(:).*y(:));
    sigmax = sqrt(sum(w(:).*(x(:)-sumx(:)).^2));
    sigmay = sqrt(sum(w(:).*(y(:)-sumy(:)).^2));
    c = sum(w(:).*(x(:)-sumx(:)).*(y(:)-sumy(:))) / (sigmax*sigmay);
end

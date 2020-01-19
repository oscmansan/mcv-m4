function disparity = stereo_computation(Ileft,Iright,minDisp,maxDisp,winSize,varargin)
%STEREO_COMPUTATION Summary of this function goes here
%   Detailed explanation goes here

cost = 'SAD';
bw = false;
asym = true;

for k=1:2:length(varargin)
    switch lower(varargin{k})
        case 'cost'
            cost = varargin{k+1};
        case 'bw'
            bw = varargin{k+1};
        case 'asym'
            asym = varargin{k+1};
        otherwise
            error(['[Unknown option ''', varargin{k}, '''.']) ;
    end
end

switch cost
    case 'SAD'
        cost = @sad;
    case 'SSD'
        cost = @ssd;
    case 'NCC'
        cost = @(x,y,w) 1-ncc(x,y,w);
    otherwise
        error('Unexpected cost type.');
end

[numRows,numCols] = size(Ileft);
r = floor(winSize/2);

gamma_col = 12;
gamma_pos = r;

disparity = zeros(size(Ileft));
tic
for i = 1+r:numRows-r
    for j = 1+r:numCols-r
        patchLeft = Ileft(i-r:i+r,j-r:j+r);

        w = ones(winSize);
        if bw
            w = bilateral_weights(patchLeft,gamma_col,gamma_pos);
        end

        minCost = Inf;
        for k = [max(1+r,j-maxDisp):max(1+r,j-minDisp) min(numCols-r,j+minDisp):min(numCols-r,j+maxDisp)]
            patchRight = Iright(i-r:i+r,k-r:k+r);

            if bw && ~asym
                w = w.*bilateral_weights(patchRight,gamma_col,gamma_pos);
            end
            w = w./sum(w(:));  % normalize weights

            c = cost(patchLeft,patchRight,w);
            if c < minCost
                minCost = c;
                idx = k;
            end
        end
        disparity(i,j) = abs(j-idx);
    end
    fprintf('%.2f %%\n',(i-r)/(numRows-2*r)*100);
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

function w = bilateral_weights(x,gammac,gammap)
    sz = size(x);
    c = ceil(sz/2);
    deltac = zeros(sz);
    deltap = zeros(sz);
    for i = 1:sz(1)
        for j = 1:sz(2)
            deltac(i,j) = abs(x(i,j)-x(c(1),c(2)));
            deltap(i,j) = sqrt((i-c(1))^2+(j-c(2))^2);
        end
    end
    w = exp(-(deltac/gammac+deltap/gammap));
end

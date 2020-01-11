function tf = iscloseall(A,B,tol,nansequal)
%tf = iscloseall(A,B) - check if two matrices are close within tolerance
%
%Usage:
%   tf = iscloseall(A,B)
%   tf = iscloseall(A,B,tol)
%   tf = iscloseall(A,B,tol,nansequal)
%
%Inputs:
%   -A: First input for comparison
%   -B: Second input for comparison
%   -tol: optional scalar tolerance to define "close".
%           {default:  max(1e4*eps(min(abs(A),abs(B))))}
%   -nansequal: logical are nans equal to themselves?  
%        Use true to replicate "iscloseall with equalnans"
%           {default: false}
%
%Outputs:
%   -tf: logical scalar, true if all values in A and B are within tol and A
%        and B are the same size and class.
%
%See also: isequal isequaln
%

    %Error checking and default assignments
    assert(any(nargin==2:4),'Two, three or four inputs expected')
    if nargin <= 3
        nansequal = false;
    else
        assert(~isempty(nansequal) && isscalar(nansequal) && ...
            islogical(nansequal) && ~isnan(nansequal),...
            'nanequal, the fourth input, is expected to be logical scalar');
    end
    if nargin == 2
        tol = 1e4*eps(min(abs(A),abs(B)));
        tol = max(tol(:));
    else
        assert(~isempty(tol) && isscalar(tol) && isnumeric(tol) && ...
            ~isnan(tol),'tol, the third input, is expected to be numeric scalar');    
    end

    %Engine:
    tf = false;                     %guilty unless proven otherwise
    if strcmp(class(A),class(B))    %classes equal
        if isequal(size(A),size(B)) %sizes first so subtraction doesn't error
            A = A(:);               %Column vectors
            B = B(:);
            ABabsdiff = abs(A-B);   %absolute differences
            if nansequal                
                %If nan in both A and B, set to 0
                ABabsdiff(isnan(A)&isnan(B)) = 0;
            end
            if all(ABabsdiff<tol)       
                %all absdiff within tol
                tf = true;
            end
        end
    end
end
function B = downsize_mean(A, factor)
% DOWNSIZE_MEAN downsizes a matrix or vector by a given factor using the mean.
%   - Works for 2D matrices (MxN), row vectors (1xN), or column vectors (Nx1).
%   - Input dimensions must be divisible by the factor.

    sz = size(A);
    
    if isvector(A)
        len = numel(A);
        if mod(len, factor) ~= 0
            error('Vector length must be divisible by the reduction factor.');
        end
        
        % Convert to row vector for consistent handling
        A_row = A(:)'; % make it 1 x N
        A_reshaped = reshape(A_row, factor, []);
        B = mean(A_reshaped, 1);
        
        % Return result in original orientation
        if iscolumn(A)
            B = B';  % convert back to column vector
        end
        
    elseif ismatrix(A)
        [m, n] = size(A);
        if mod(m, factor) ~= 0 || mod(n, factor) ~= 0
            error('Matrix dimensions must be divisible by the reduction factor.');
        end
        
        B = mean(reshape(A, factor, m/factor, factor, n/factor), [1 3]);
        B = squeeze(B);
        
    else
        error('Input must be a 1D or 2D matrix.');
    end
end

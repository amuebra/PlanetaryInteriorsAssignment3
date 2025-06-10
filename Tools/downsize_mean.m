function B = downsize_mean(A, factor)
% downsize_mean downsizes a matrix by a given factor. It uses the
% arthmetic mean to average the multiple input cells. The script also checks if
% the input dimensions are divisible by the given factor, before the
% downsizing takes place.

    % Check input validity
    [m, n] = size(A);
    if mod(m, factor) ~= 0 || mod(n, factor) ~= 0
        error('Matrix dimensions must be divisible by the reduction factor.');
    end

    % Reshape and average
    B = mean(reshape(A, factor, m/factor, factor, n/factor), [1 3]);
    B = squeeze(B);
end
function X = B_XV(V, m)
    [n, r] = size(V);
    B_blocks = cell(1, r);

    for j = 1:r
        v_j = V(:, j);   % n\times 1
        B_blocks{j} = kron(eye(m), v_j);  %(m*n)\times m
    end

    X = blkdiag(B_blocks{:});
end
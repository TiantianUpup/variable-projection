function X = B_UY(U, n)
    [m, r] = size(U);
    B_blocks = cell(1, r);

    for j = 1:r
        u_j = U(:, j);  % m\times 1
        B_blocks{j} = kron(u_j, eye(n)); % (m*n)\times m
    end

    X = blkdiag(B_blocks{:});
end
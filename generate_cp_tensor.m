function X = generate_cp_tensor(U, V, W)
    [m,r] = size(U);
    [n,r] = size(V);
    [p,r] = size(W);
    
    X = zeros(m, n, p);
    
    for i = 1:r
        u = U(:, i);
        v = V(:, i);
        w = W(:, i);
        
        rank1_tensor = zeros(m, n, p);
        for a = 1:m
            for b = 1:n
                for c = 1:p
                    rank1_tensor(a, b, c) = u(a) * v(b) * w(c);
                end
            end
        end
        
        X = X + rank1_tensor;
    end
end
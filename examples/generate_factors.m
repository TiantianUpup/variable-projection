function [U, V, W] = generate_factors(m, n, p, R)
    U = gen_factor(m, R);
    V = gen_factor(n, R);
    W = gen_factor(p, R);
end

function A = gen_factor(dim, R)
    A = zeros(dim, R);
    for r = 1:R
        k = round(dim * 0.1);
        idx = randperm(dim, k);
        col = zeros(dim, 1);
        col(idx) = 100 * rand(k, 1);
        col(setdiff(1:dim, idx)) = rand(dim - k, 1);
        A(:, r) = col / sum(col);
    end
end
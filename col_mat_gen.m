function X = col_mat_gen(dim,r,c)
    % Step 1: generate an r\times r matrix K that has diagonal elements 1 and off-diagonal elements c
    K = ones(r)*c + (1-c)*eye(r);
    C = chol(K);
    
    X=cell(3);
    % Step 2

    for i=1:3
        T = rand(dim,r);
        [Q, ~] = qr(T,0);

        X{i} = Q*C;
    end
end   
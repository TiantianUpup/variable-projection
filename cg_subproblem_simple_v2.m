function runhist = cg_subproblem_simple_v2(Aparas,cg_paras,paras,B,gamma)
    timect = cputime;

    % convergence condition
    itmax=cg_paras.itmax;
    tol=cg_paras.tol;
    
    % initialization
    r=paras.r;
    m=paras.m;
    n=paras.n;
    X=zeros(m+n, r); % direction D=[A;B];
    %X=X0;
    % calculation of R
    R = B;

    iter = 0;
    R_norm = norm(R(:));

    while (iter < itmax && R_norm > tol)
        % step 1
        iter = iter + 1;

        % step 2
        if iter == 1
            C = R;
        else
            beta = norm(R(:))^2 / norm(R_pre(:))^2;
            C = R + beta * C;
        end

        % Step 3
        alpha = norm (R(:))^2 / trace(C' * Aope_simple_v2(Aparas, paras, C,gamma));

        % Step 4
        X = X + alpha * C;

        % Step 5
        R_pre = R;
        R = R - alpha * Aope_simple_v2(Aparas, paras, C, gamma);

        R_norm = norm(R(:));
    end

    cput = cputime - timect;
    %fprintf("cg method iter=%d, R_norm=%3.7f, time cost is %3.4f\n",iter,norm(R_norm(:)),cput);
    runhist.R_norm = R_norm;
    runhist.X = X;
    runhist.iter = iter;
end

function runhist = cg_hy_subproblem (Aparas,cg_paras,paras,B, gamma)
    timect = cputime;

    % convergence condition
    itmax=cg_paras.itmax;
    tol=cg_paras.tol;
    mu=cg_paras.mu;
    lambda=cg_paras.lambda;

    % initialization
    r=paras.r;
    m=paras.m;
    n=paras.n;
    x=zeros((m+n)*r,1); % direction D=[A;B];
    %X=X0;
    % calculation of R
    r = B(:);

    iter = 0;
    r_norm = norm(r);

    while (iter < itmax && r_norm > tol)
        % step 1
        iter = iter + 1;

        % step 2
        if iter == 1
            c = r;
        else
            beta = norm(r)^2 / norm(r_pre)^2;
            c = r + beta * c;
        end

        % Step 3
        alpha = norm (r)^2 / (c' * Aope_hy(Aparas, paras, c, gamma,mu,lambda));

        % Step 4
        x = x + alpha * c;

        % Step 5
        r_pre = r;
        r = r - alpha * Aope_hy(Aparas, paras, c, gamma,mu,lambda);

        r_norm = norm(r);
    end

    cput = cputime - timect;
    fprintf("cg method iter=%d, r_norm=%3.7f, time cost is %3.4f\n",iter,norm(r_norm(:)),cput);
    runhist.r_norm = r_norm;
    runhist.x = x;
    runhist.iter = iter;
end

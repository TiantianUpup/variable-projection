function runhist = cg_subproblem (Aparas, paras, B, gamma,lambda)
    %fprintf("lambda is %3.4f\n",lambda);
    % convergence condition
    %itmax=paras.itmax;
    itmax = 500;
    %tol=paras.tol;
    tol = 1e-6;

    % initialization
    r = paras.r;
    m=paras.m;
    n=paras.n;
    X = zeros(m+n, r);
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
        alpha = norm (R(:))^2 / trace(C' * Aope(Aparas, paras, C, gamma,lambda));

        % Step 4
        X = X + alpha * C;

        % Step 5
        R_pre = R;
        R = R - alpha * Aope(Aparas, paras, C, gamma,lambda);

        R_norm = norm(R(:));
    end

    fprintf("cg method iter=%d, R_norm=%3.7f\n",iter,norm(R_norm(:)));
    runhist.R_norm = R_norm;
    runhist.X = X;
    runhist.iter = iter;
end

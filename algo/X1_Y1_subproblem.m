function runhist = X1_Y1_subproblem (Aparas, B, cg_paras, paras)
    % initialization
    r=paras.r;
    X = zeros(2*r,r);
    R=B;

    iter = 0;
    R_norm = norm(R(:));
    itmax=cg_paras.itmax;
    tol=cg_paras.tol;

    while (iter  < itmax && R_norm>tol)
        % step 1
        iter = iter + 1;

        % step 2
        if iter == 1
            C = R;
        else
            beta = norm(R(:))^2/norm(R_pre(:))^2;
            C=R+beta*C;
        end

        % Step 3
        alpha = norm (R(:))^2 / trace(C'*AadjointA(Aparas, paras, C));

        % Step 4
        X = X + alpha * C;

        % Step 5
        R_pre=R;
        R = R - alpha*AadjointA(Aparas, paras, C);

        R_norm=norm(R(:));
    end

    runhist.R_norm=R_norm;
    runhist.X=X;
    runhist.iter=iter;
end    
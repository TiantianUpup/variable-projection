function runhist = X1_Y1_subproblem (Aparas, B, paras)
    gamma=0;

    % initialization
    r=paras.r;
    X = zeros(2*r,r);
    R=B;

    iter = 0;
    R_norm = norm(R(:));
    itmax=500;
    tol=1e-6;
    %fprintf("X1_Y1_subproblem R norm is %3.4f\n", R_norm);

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
        alpha = norm (R(:))^2 / trace(C'*AadjointA(Aparas, paras, C,gamma));

        % Step 4
        X = X + alpha * C;

        % Step 5
        R_pre=R;
        R = R - alpha*AadjointA(Aparas, paras, C,gamma);

        R_norm=norm(R(:));
    end

    %fprintf("X1_Y1_subproblem iter=%d\n",iter);

    runhist.R_norm=R_norm;
    runhist.X=X;
    runhist.iter=iter;
    %fprintf("------------------ X_1_Y1_subproblem iter=%d ---------------------------\n",iter);
end    